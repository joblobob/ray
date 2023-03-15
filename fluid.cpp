#include "fluid.h"

#include <QGraphicsView>
#include <QOpenGLWidget>
#include <QPainter>
#include <QTimer>

#include <execution>

#include "ui_fluid.h"


import Obstacle;
import Constants;
import format;

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
	m_f = FlipFluid(constants::maxwidth, constants::maxheight, constants::cellHeight, constants::particleRadius, constants::maxParticles);

	// create particles
	int count = 0;
	for (int i = 0; i < constants::numX; i++) {
		for (int j = 0; j < constants::numY; j++, count++) {
			const float posX       = constants::cellHeight + constants::particleRadius + constants::dx * i + (j % 2 == 0 ? 0.0f : constants::particleRadius);
			const float posY       = constants::cellHeight + constants::particleRadius + constants::dy * j;
			m_f.particleMap[count] = { count, posX, posY };
		}
	}


	// setup grid cells for tank

	float n = constants::fNumY;
	count   = 0;
	for (int i = 0; i < constants::fNumX; i++) {
		for (int j = 0; j < constants::fNumY; j++) {
			float s = 1.0f; // fluid
			if (i == 0 || i == constants::fNumX - 1 || j == 0)
				s = 0.0f; // solid
			m_f.gridCells[i * n + j].s        = s;
			m_f.gridCells[i * n + j].cellNumX = i;
			m_f.gridCells[i * n + j].cellNumY = j;
			count++;
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

//draw to pixmap
void fluid::drawimage()
{
	if (constants::showParticles) {
		m_pointsPixmap.fill(Qt::black);

		float pointSize       = 2.0f * constants::particleRadius;
		auto setParticleColor = [&](auto&& item) {
			m_painter->setBrush({ QColor::fromRgbF(item.colorR, item.colorG, item.colorB) });
			m_painter->setPen({ QColor::fromRgbF(item.colorR, item.colorG, item.colorB) });
			m_painter->drawEllipse(item.posX, item.posY, pointSize, pointSize);
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