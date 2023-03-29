#include "fluid.h"

#include <QGraphicsView>
#include <QOpenGLWidget>
#include <QPainter>
#include <QTimer>

#include <execution>

#include <mdspan.hpp>

#include "ui_fluid.h"


import Obstacle;
import Constants;
import format;
//import mdspan;


fluid::fluid(QWidget* parent) :
    ui(new Ui::fluid), m_paused(true), m_pointsPixmap(constants::maxwidth, constants::maxheight), m_pixmapSceneItem(nullptr)
{
	ui->setupUi(this);

	m_pointsPixmap.fill(Qt::black);
	m_pixmapScene = new QGraphicsScene(this);
	m_painter     = new QPainter(&m_pointsPixmap);

	QOpenGLWidget* gl = new QOpenGLWidget();
	ui->m_graphView->setViewport(gl);
	ui->m_graphView->setRenderHint(QPainter::Antialiasing, false);
	ui->m_graphView->setOptimizationFlag(QGraphicsView::DontSavePainterState, true);
	ui->m_graphView->setOptimizationFlag(QGraphicsView::DontAdjustForAntialiasing, true);
	ui->m_graphView->setOptimizationFlag(QGraphicsView::IndirectPainting, true);
	ui->m_graphView->setViewportUpdateMode(QGraphicsView::FullViewportUpdate);
	ui->m_graphView->setCacheMode(QGraphicsView::CacheBackground);

	ui->m_graphView->setScene(m_pixmapScene);
	ui->m_graphView->scale(1.0f, -1.0f); // invert Y axis
	ui->m_graphView->setSceneRect(0.0f, 0.0f, constants::maxwidth, constants::maxheight);


	QTimer::singleShot(0, this, &fluid::update);
}

void fluid::setupScene()
{
	// create fluid
	m_f = FlipFluid(constants::maxwidth, constants::maxheight, constants::cellHeight, constants::maxParticles);

	std::experimental::mdspan particleView(m_f.particleMap.data(), constants::numX, constants::numY);
	// create particles
	int count = 0;
	for (std::size_t i = 0; i < particleView.extent(0); i++) {
		for (std::size_t j = 0; j < particleView.extent(1); j++, count++) {
			float rad = constants::particleDefaultRadius;
			if (count < 3000)
				rad *= 4.0f;
			const float posX   = constants::cellHeight + rad + constants::dx * i + (j % 2 == 0 ? 0.0f : rad);
			const float posY   = constants::cellHeight + rad + constants::dy * j;
			particleView(i, j) = { .id = count, .posX = posX, .posY = posY, .radius = rad };
		}
	}



	// setup grid cells for tank
	std::experimental::mdspan gridView(m_f.gridCells.data(), constants::fNumX, constants::fNumY);
	float n = constants::fNumY;
	for (std::size_t i = 0; i < gridView.extent(0); i++) {
		for (std::size_t j = 0; j < gridView.extent(1); j++) {
			float s = 1.0f; // fluid
			if (i == 0 || i == constants::fNumX - 1 || j == 0)
				s = 0.0f; // solid
			gridView(i, j).s        = s;
			gridView(i, j).cellNumX = i;
			gridView(i, j).cellNumY = j;
		}
	}


	setupObstacle(m_f.gridCells, 0.0f, 0.0f, m_f.obstacle, true);
	m_obstacleItem =
	    m_pixmapScene->addEllipse(0, 0, ObstacleConstants::radius * 2.0f, ObstacleConstants::radius * 2.0f, QPen(Qt::NoPen), QBrush(Qt::red));
	m_obstacleItem->setFlag(QGraphicsItem::ItemIsMovable);
	m_obstacleItem->setFlag(QGraphicsItem::ItemIgnoresParentOpacity, true);
}

void fluid::simulate()
{
	if (!m_paused) {
		auto result = m_f.simulate(constants::dt, true);
		if (result.size() > 0) {
			constexpr std::string msgFormat { "{}: {:.2f} us" };
			ui->label_2->setText(std::format(msgFormat, result[0].message, result[0].nsElapsed * 0.001f).data());
			ui->label_3->setText(std::format(msgFormat, result[1].message, result[1].nsElapsed * 0.001f).data());
			ui->label_4->setText(std::format(msgFormat, result[2].message, result[2].nsElapsed * 0.001f).data());
			ui->label_5->setText(std::format(msgFormat, result[3].message, result[3].nsElapsed * 0.001f).data());
			ui->label_6->setText(std::format(msgFormat, result[4].message, result[4].nsElapsed * 0.001f).data());
			ui->label_7->setText(std::format(msgFormat, result[5].message, result[5].nsElapsed * 0.001f).data());
			ui->label_8->setText(std::format(msgFormat, result[6].message, result[6].nsElapsed * 0.001f).data());
			ui->label_9->setText(std::format(msgFormat, result[7].message, result[7].nsElapsed * 0.001f).data());
			ui->label_10->setText(std::format(msgFormat, result[8].message, result[8].nsElapsed * 0.001f).data());
		}
	}
}


QBrush soilBrush { QColor { 161, 103, 23 } };
QPen soilPen { QColor { 161, 103, 23 } };
QBrush lightSoilBrush { QColor { 191, 133, 53 } };
QPen lightSoilPen { QColor { 191, 133, 53 } };
QBrush rockBrush { Qt::darkGray };
QPen rockPen { Qt::darkGray };
QBrush grassBrush { Qt::darkGreen };
QPen grassPen { Qt::darkGreen };

//draw to pixmap
void fluid::drawimage()
{
	if (constants::showParticles) {
		m_pointsPixmap.fill(Qt::black);

		auto setParticleColor = [&](const Particle& particle) {
			auto color = QColor::fromRgbF(particle.colorR, particle.colorG, particle.colorB);
			if (particle.id < 3000) {
				m_painter->setBrush(rockBrush);
				m_painter->setPen(rockPen);
			} else if (particle.id < 10000) {
				m_painter->setBrush(soilBrush);
				m_painter->setPen(soilPen);
			} else if (particle.id < 15000) {
				m_painter->setBrush(lightSoilBrush);
				m_painter->setPen(lightSoilPen);
			} else if (particle.id < 17000) {
				m_painter->setBrush(grassBrush);
				m_painter->setPen(grassPen);
			} else {
				m_painter->setBrush({ color });
				m_painter->setPen({ color });
			}
			m_painter->drawEllipse(particle.posX, particle.posY, particle.radius * 2.0f, particle.radius * 2.0f);
		};

		std::ranges::for_each(m_f.particleMap, setParticleColor);
	}

	//obstacle
	m_painter->setBrush(Qt::red);
	m_painter->setPen(Qt::red);
	m_painter->drawEllipse(
	    ((float)m_obstacleItem->x()), ((float)m_obstacleItem->y()), ObstacleConstants::radius * 2.0f, ObstacleConstants::radius * 2.0f);


	drawImageToScene();
}

void fluid::update()
{
	float totalElapsed { m_timer.nsecsElapsed() * 0.000001f };


	m_simtimer.restart();

	setupObstacle(m_f.gridCells,
	    ((float)m_obstacleItem->x() + ObstacleConstants::radius),
	    ((float)m_obstacleItem->y() + ObstacleConstants::radius),
	    m_f.obstacle,
	    false);

	simulate();
	float simElapsed { m_simtimer.nsecsElapsed() * 0.000001f };
	m_simtimer.restart();
	//draw();
	drawimage();

	float simDrawElapsed { m_simtimer.nsecsElapsed() * 0.000001f };
	fps++;

	std::string fpsLabel = std::format("Frame: {} Simulate: {:.2f} ms SimDraw: {:.2f} ms roundtrip: {:.2f} ms total {:.2f}",
	    fps,
	    simElapsed,
	    simDrawElapsed,
	    totalElapsed,
	    simElapsed + simDrawElapsed + totalElapsed);
	ui->label->setText(fpsLabel.data());

	m_timer.restart();

	QTimer::singleShot(0, this, &fluid::update);
}


void fluid::drawImageToScene()
{
	if (m_pixmapSceneItem == nullptr) {
		m_pixmapSceneItem = m_pixmapScene->addPixmap(m_pointsPixmap);
		m_pixmapScene->setSceneRect(m_pointsPixmap.rect());
	} else
		m_pixmapSceneItem->setPixmap(m_pointsPixmap);
}