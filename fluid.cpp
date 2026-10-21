#include "fluid.h"

#include <QGraphicsView>
#ifndef __EMSCRIPTEN__
#include <QOpenGLWidget>
#endif
#include <QPainter>
#include <QTimer>

#include <execution>
#include <format>

#include <mdspan>
#include <vector>

#include "ui_fluid.h"


import Obstacle;
import Constants;
import FlipFluid;
import BaseStructures;
//import mdspan;

struct fluid::Impl {
	FlipFluid f;
};

class ParticleOverlayItem final : public QGraphicsItem
{
public:
	explicit ParticleOverlayItem(const std::vector<Particle>* particles) : m_particles(particles)
	{
		setZValue(1.0);
	}

	QRectF boundingRect() const override
	{
		return { 0.0, 0.0, constants::maxwidth, constants::maxheight };
	}

	void paint(QPainter* painter, const QStyleOptionGraphicsItem*, QWidget*) override
	{
		if (!constants::showParticles || m_particles == nullptr)
			return;

		// Paint in viewport pixels after using the scene transform solely to
		// map each simulated position. This keeps particle diameter independent
		// of fitInView's zoom level.
		const QTransform sceneToViewport = painter->worldTransform();
		painter->save();
		painter->resetTransform();
		painter->setPen(Qt::NoPen);

		constexpr qreal diameter = constants::dx;
		constexpr qreal radius   = diameter * 0.5;
		for (const Particle& item : *m_particles) {
			QColor color;
			if (item.id < 1500)
				color = Qt::darkGray;
			else if (item.id < 5000)
				color = { 161, 103, 23 };
			else if (item.id < 7500)
				color = { 191, 133, 53 };
			else if (item.id < 8000)
				color = Qt::darkGreen;
			else
				color = QColor::fromRgbF(item.colorR, item.colorG, item.colorB);

			const QPointF position = sceneToViewport.map(QPointF(item.posX, item.posY));
			painter->setBrush(color);
			painter->drawEllipse(position.x() - radius, position.y() - radius, diameter, diameter);
		}
		painter->restore();
	}

private:
	const std::vector<Particle>* m_particles;
};

fluid::fluid(QWidget* parent) :
    ui(new Ui::fluid), m_paused(true), m_impl(std::make_unique<Impl>()), m_pointsPixmap(constants::maxwidth, constants::maxheight), m_pixmapSceneItem(nullptr)
{
	ui->setupUi(this);

	m_pointsPixmap.fill(Qt::black);
	m_pixmapScene = new QGraphicsScene(this);

#ifndef __EMSCRIPTEN__
	QOpenGLWidget* gl = new QOpenGLWidget();
	ui->m_graphView->setViewport(gl);
#endif
	ui->m_graphView->setRenderHint(QPainter::Antialiasing, false);
	ui->m_graphView->setOptimizationFlag(QGraphicsView::DontSavePainterState, true);
	ui->m_graphView->setOptimizationFlag(QGraphicsView::DontAdjustForAntialiasing, true);
	ui->m_graphView->setOptimizationFlag(QGraphicsView::IndirectPainting, true);
	ui->m_graphView->setViewportUpdateMode(QGraphicsView::FullViewportUpdate);
	ui->m_graphView->setCacheMode(QGraphicsView::CacheBackground);

	ui->m_graphView->setScene(m_pixmapScene);
	ui->m_graphView->scale(1.0f, -1.0f); // invert Y axis
	ui->m_graphView->setSceneRect(0.0f, 0.0f, constants::maxwidth, constants::maxheight);


	// The simulation advances by a fixed 1/60 s per update. Pace updates to
	// that rate so simulation speed does not depend on how fast this loop runs.
	QTimer::singleShot(16, this, &fluid::update);
}

fluid::~fluid()
{
	delete ui;
}

void fluid::setupScene()
{
	// create fluid
	m_impl->f = FlipFluid(constants::maxwidth, constants::maxheight, constants::cellHeight, constants::maxParticles);

	std::mdspan particleView(m_impl->f.particleMap.data(), constants::numX, constants::numY);
	// create particles
	int count = 0;
	for (std::size_t i = 0; i < particleView.extent(0); i++) {
		for (std::size_t j = 0; j < particleView.extent(1); j++, count++) {
			const float posX   = constants::cellHeight + constants::particleRadius + constants::dx * i + (j % 2 == 0 ? 0.0f : constants::particleRadius);
			const float posY   = constants::cellHeight + constants::particleRadius + constants::dy * j;
			particleView[i, j] = { count, posX, posY };
		}
	}


	// setup grid cells for tank
	std::mdspan gridView(m_impl->f.gridCells.data(), constants::fNumX, constants::fNumY);
	float n = constants::fNumY;
	for (std::size_t i = 0; i < gridView.extent(0); i++) {
		for (std::size_t j = 0; j < gridView.extent(1); j++) {
			float s = 1.0f; // fluid
			if (i == 0 || i == constants::fNumX - 1 || j == 0)
				s = 0.0f; // solid
			gridView[i, j].s        = s;
			gridView[i, j].cellNumX = i;
			gridView[i, j].cellNumY = j;
		}
	}


	setupObstacle(m_impl->f.gridCells, 0.0f, 0.0f, m_impl->f.obstacle, false);
	drawImageToScene();
	m_pixmapSceneItem->setZValue(-1.0);

	if (m_obstacleItem != nullptr)
		delete m_obstacleItem;
	m_obstacleItem =
	    m_pixmapScene->addEllipse(0, 0, ObstacleConstants::radius * 2.0f, ObstacleConstants::radius * 2.0f, QPen(Qt::NoPen), QBrush(Qt::red));
	m_obstacleItem->setFlag(QGraphicsItem::ItemIsMovable);
	m_obstacleItem->setFlag(QGraphicsItem::ItemIgnoresParentOpacity, true);
	m_obstacleItem->setZValue(2.0);

	if (m_particleOverlayItem == nullptr) {
		m_particleOverlayItem = new ParticleOverlayItem(&m_impl->f.particleMap);
		m_pixmapScene->addItem(m_particleOverlayItem);
	}
	m_pixmapScene->setSceneRect(m_pointsPixmap.rect());
	ui->m_graphView->fitInView(m_pixmapScene->sceneRect(), Qt::KeepAspectRatio);
	// setupScene runs before the dialog is shown, when the viewport can still
	// have its designer-time size. Refit after Qt has completed its layout.
	QTimer::singleShot(0, this, [this] {
		ui->m_graphView->fitInView(m_pixmapScene->sceneRect(), Qt::KeepAspectRatio);
		if (m_particleOverlayItem != nullptr)
			m_particleOverlayItem->update();
	});
}

void fluid::resizeEvent(QResizeEvent* event)
{
	QDialog::resizeEvent(event);
	if (m_pixmapScene != nullptr)
		ui->m_graphView->fitInView(m_pixmapScene->sceneRect(), Qt::KeepAspectRatio);
}

void fluid::simulate()
{
	if (!m_paused) {
		auto result = m_impl->f.simulate(constants::dt, false);
		if (result.size() > 0) {
			//constexpr static std::string msgFormat { "{}: {:.2f} us" };
			ui->label_2->setText(std::format("{}: {:.2f} us", result[0].message, result[0].nsElapsed * 0.001f).data());
			ui->label_3->setText(std::format("{}: {:.2f} us", result[1].message, result[1].nsElapsed * 0.001f).data());
			ui->label_4->setText(std::format("{}: {:.2f} us", result[2].message, result[2].nsElapsed * 0.001f).data());
			ui->label_5->setText(std::format("{}: {:.2f} us", result[3].message, result[3].nsElapsed * 0.001f).data());
			ui->label_6->setText(std::format("{}: {:.2f} us", result[4].message, result[4].nsElapsed * 0.001f).data());
			ui->label_7->setText(std::format("{}: {:.2f} us", result[5].message, result[5].nsElapsed * 0.001f).data());
			ui->label_8->setText(std::format("{}: {:.2f} us", result[6].message, result[6].nsElapsed * 0.001f).data());
			ui->label_9->setText(std::format("{}: {:.2f} us", result[7].message, result[7].nsElapsed * 0.001f).data());
			ui->label_10->setText(std::format("{}: {:.2f} us", result[8].message, result[8].nsElapsed * 0.001f).data());
		}
	}
}


void fluid::drawimage()
{
	if (m_particleOverlayItem != nullptr)
		m_particleOverlayItem->update();
}

void fluid::update()
{
	float totalElapsed { m_timer.nsecsElapsed() * 0.000001f };


	m_simtimer.restart();

	setupObstacle(m_impl->f.gridCells,
	    ((float)m_obstacleItem->x() + ObstacleConstants::radius),
	    ((float)m_obstacleItem->y() + ObstacleConstants::radius),
	    m_impl->f.obstacle,
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

	QTimer::singleShot(16, this, &fluid::update);
}


void fluid::drawImageToScene()
{
	if (m_pixmapSceneItem == nullptr) {
		m_pixmapSceneItem = m_pixmapScene->addPixmap(m_pointsPixmap);
		m_pixmapSceneItem->setCacheMode(QGraphicsItem::DeviceCoordinateCache);
		m_pixmapScene->setSceneRect(m_pointsPixmap.rect());
	}
}
