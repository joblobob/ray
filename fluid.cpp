#include "fluid.h"

#include <QGraphicsView>
#include <QTimer>
#include <execution>

#include "ui_fluid.h"


import Obstacle;
import Constants;
import format;

fluid::fluid(QWidget* parent) :
    ui(new Ui::fluid), m_paused(true), m_imageCanvas(constants::maxwidth, constants::maxheight, QImage::Format_RGB32), m_pixmapSceneItem(nullptr)
{
	ui->setupUi(this);

	m_imageCanvas.fill(Qt::black);
	m_pixmapScene = new QGraphicsScene(this);

	m_scene = new QGraphicsScene(this);

	ui->m_graphView->setRenderHint(QPainter::Antialiasing, false);
	ui->m_graphView->setOptimizationFlag(QGraphicsView::DontSavePainterState, true);
	ui->m_graphView->setOptimizationFlag(QGraphicsView::DontAdjustForAntialiasing, true);
	ui->m_graphView->setOptimizationFlag(QGraphicsView::IndirectPainting, true);
	ui->m_graphView->setViewportUpdateMode(QGraphicsView::FullViewportUpdate);
	ui->m_graphView->setCacheMode(QGraphicsView::CacheBackground);

	m_scene->setItemIndexMethod(QGraphicsScene::NoIndex);

	ui->m_graphView->setScene(m_scene);
	ui->m_graphView->scale(1.0f, -1.0f); // invert Y axis
	ui->m_graphView->setSceneRect(0.0f, 0.0f, constants::maxwidth, constants::maxheight);

	m_drawTimer.setInterval(16);
	connect(&m_drawTimer, &QTimer::timeout, this, &fluid::update);
	m_drawTimer.start();
}

void fluid::setupScene()
{
	m_scene->clear();
	m_particleItems.clear();
	m_gridItems.clear();


	// create fluid
	m_f = FlipFluid(constants::maxwidth, constants::maxheight, constants::cellHeight, constants::particleRadius, constants::maxParticles);

	// create particles
	float r2 = constants::particleRadius * 2.0f;

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

	float sizeRect = constants::cellHeight;

	//setup grid
	int nbGridItems = constants::fNumX * constants::fNumY;
	m_gridItems.reserve(nbGridItems);

	auto pointSizeGrid = constants::cellHeight;
	if (constants::showGrid) {
		auto gridPen         = QPen(QColor(127, 127, 127));
		int yesMyIndexAtHand = 0;

		auto setGridItems = [&](const Cell& cell) {
			auto item = m_scene->addRect(0, 0, pointSizeGrid, pointSizeGrid, gridPen, QBrush(QColor::fromRgbF(cell.colorR, cell.colorG, cell.colorB)));
			m_gridItems.push_back(item);
			item->setPos(cell.cellNumX * pointSizeGrid, cell.cellNumY * pointSizeGrid);
			item->setData(0, yesMyIndexAtHand);
			item->setFlag(QGraphicsItem::ItemIgnoresParentOpacity, true);
			yesMyIndexAtHand++;
		};
		std::for_each(std::execution::seq, m_f.gridCells.begin(), m_f.gridCells.end(), setGridItems);
	}

	//setup particles of water
	m_particleItems.reserve(constants::maxParticles);
	if (constants::showParticles) {
		float pointSize = 2.0f * constants::particleRadius;
		for (int i = 0; i < constants::maxParticles; i++) {}
		auto setParticleItems = [&](const Particle& particle) {
			auto particleColor = QColor::fromRgbF(particle.colorR, particle.colorG, particle.colorB);
			m_particleItems.push_back(m_scene->addEllipse(0.0f, 0.0f, pointSize, pointSize, QPen(Qt::NoPen), QBrush(particleColor)));

			m_particleItems[particle.id]->setFlag(QGraphicsItem::ItemIgnoresParentOpacity, true);
			m_particleItems[particle.id]->setPos(particle.posX - constants::particleRadius, particle.posY - constants::particleRadius);
			m_particleItems[particle.id]->setData(0, particle.id);
		};
		std::for_each(std::execution::seq, m_f.particleMap.cbegin(), m_f.particleMap.cend(), setParticleItems);
	}


	setupObstacle(m_f.gridCells, 0.0f, 0.0f, m_f.obstacle, true);
	m_obstacleItem = m_scene->addEllipse(0, 0, ObstacleConstants::radius * 2.0f, ObstacleConstants::radius * 2.0f, QPen(Qt::NoPen), QBrush(Qt::red));
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

void fluid::draw()
{
	//grid
	if (constants::showGrid) {
		auto setCellColor = [&](QGraphicsRectItem* item) {
			auto i          = item->data(0).toInt();
			const Cell cell = m_f.gridCells[i];
			item->setBrush(QBrush(QColor::fromRgbF(cell.colorR, cell.colorG, cell.colorB)));
		};

		std::for_each(std::execution::par_unseq, m_gridItems.begin(), m_gridItems.end(), setCellColor);
	}

	//particles
	//setup particles of water
	if (constants::showParticles) {
		auto setParticleColor = [&](QGraphicsEllipseItem* item) {
			auto i             = item->data(0).toInt();
			auto particleColor = QColor::fromRgbF(m_f.particleMap[i].colorR, m_f.particleMap[i].colorG, m_f.particleMap[i].colorB);
			item->setPos(m_f.particleMap[i].posX - constants::particleRadius, m_f.particleMap[i].posY - constants::particleRadius);
			item->setBrush(QBrush(particleColor));
		};

		std::for_each(std::execution::par_unseq, m_particleItems.begin(), m_particleItems.end(), setParticleColor);
	}
}



void fluid::update()
{
	std::string fpsLabel = std::format("Simulate: {} fps Draw: {:.2f} ms", fps, m_timer.nsecsElapsed() * 0.000001f);
	ui->label->setText(fpsLabel.data());

	m_timer.restart();

	setupObstacle(m_f.gridCells,
	    ((float)m_obstacleItem->x() + ObstacleConstants::radius),
	    ((float)m_obstacleItem->y() + ObstacleConstants::radius),
	    m_f.obstacle,
	    false);

	simulate();
	draw();
	fps++;
}

void fluid::drawImageToScene()
{
	if (m_pixmapSceneItem == nullptr) {
		m_pixmapSceneItem = m_pixmapScene->addPixmap(QPixmap::fromImage(m_imageCanvas));
		m_pixmapScene->setSceneRect(m_imageCanvas.rect());
	} else
		m_pixmapSceneItem->setPixmap(QPixmap::fromImage(m_imageCanvas));
}