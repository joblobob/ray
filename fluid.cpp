#include "fluid.h"

#include <QGraphicsView>
#ifndef __EMSCRIPTEN__
#include <QOpenGLWidget>
#endif
#include <QCoreApplication>
#include <QDateTime>
#include <QFile>
#include <QHash>
#include <QPainter>
#include <QPainterPath>
#include <QTextStream>
#include <QTimer>

#include <array>
#include <format>
#include <utility>

#include <mdspan>
#include <vector>

#include "ui_fluid.h"


import Obstacle;
import Constants;
import FlipFluid;
import BaseStructures;
//import mdspan;

namespace {
constexpr bool performanceLoggingEnabled = false;

#ifndef __EMSCRIPTEN__
class PerformanceLog final
{
public:
	PerformanceLog() : m_file(QCoreApplication::applicationDirPath() + "/fluid_perf.log"), m_stream(&m_file)
	{
		if (m_file.open(QIODevice::WriteOnly | QIODevice::Truncate | QIODevice::Text)) {
			m_stream << "utc_time,event,details\n";
			m_stream.flush();
		}
	}

	~PerformanceLog()
	{
		m_stream.flush();
		m_file.flush();
	}

	void append(const QString& event, const QString& details)
	{
		if (!m_file.isOpen())
			return;

		const QString timestamp = QDateTime::currentDateTimeUtc().toString(Qt::ISODateWithMs);
		m_stream << timestamp << ',' << event << ',' << details << '\n';
		if (++m_eventsSinceFlush >= 32) {
			m_stream.flush();
			m_file.flush();
			m_eventsSinceFlush = 0;
		}
	}

private:
	QFile m_file;
	QTextStream m_stream;
	unsigned int m_eventsSinceFlush = 0;
};

PerformanceLog& performanceLog()
{
	static PerformanceLog log;
	return log;
}

void resetPerformanceLog()
{
	if constexpr (performanceLoggingEnabled)
		(void)performanceLog();
}

void logPerformanceEvent(const QString& event, const QString& details)
{
	Q_UNUSED(event);
	Q_UNUSED(details);
	if constexpr (performanceLoggingEnabled)
		performanceLog().append(event, details);
}
#else
void resetPerformanceLog() {}

void logPerformanceEvent(const QString& event, const QString& details)
{
	Q_UNUSED(event);
	Q_UNUSED(details);
	if constexpr (performanceLoggingEnabled) {
		const QString timestamp = QDateTime::currentDateTimeUtc().toString(Qt::ISODateWithMs);
		qInfo().noquote() << timestamp << event << details;
	}
}
#endif
} // namespace

struct fluid::Impl {
	FlipFluid f;
};

class ParticleOverlayItem final : public QGraphicsItem
{
public:
	explicit ParticleOverlayItem(const std::vector<Particle>* particles) : m_particles(particles)
	{
		setZValue(1.0);
		for (QPainterPath& path : m_fixedColorPaths)
			path.setFillRule(Qt::WindingFill);
		m_dynamicColorIndices.reserve(256);
		m_dynamicColorPaths.reserve(256);
	}

	QRectF boundingRect() const override
	{
		return { 0.0, 0.0, constants::maxwidth, constants::maxheight };
	}

	void paint(QPainter* painter, const QStyleOptionGraphicsItem*, QWidget*) override
	{
		if (!constants::showParticles || m_particles == nullptr)
			return;
		QElapsedTimer paintTimer;
		if constexpr (performanceLoggingEnabled)
			paintTimer.start();

		// Paint in viewport pixels after using the scene transform solely to
		// map each simulated position. This keeps particle diameter independent
		// of fitInView's zoom level.
		const QTransform sceneToViewport = painter->worldTransform();
		painter->save();
		painter->resetTransform();
		painter->setPen(Qt::NoPen);
		for (QPainterPath& path : m_fixedColorPaths)
			path.clear();
		for (auto& [color, path] : m_dynamicColorPaths) {
			Q_UNUSED(color);
			path.clear();
		}
		m_activeDynamicColorPaths = 0;
		if (++m_colorGeneration == 0 || m_dynamicColorIndices.size() > 4096) {
			m_dynamicColorIndices.clear();
			m_dynamicColorIndices.reserve(256);
			m_colorGeneration = 1;
		}

		constexpr qreal diameter = constants::dx;
		constexpr qreal radius   = diameter * 0.5;
		const std::array<QColor, 4> fixedColors {
		    QColor(Qt::darkGray), QColor(161, 103, 23), QColor(191, 133, 53), QColor(Qt::darkGreen)
		};

		// Build all batches in one particle pass. The fixed-color IDs precede
		// dynamic-color IDs, so they are still painted first below.
		for (const Particle& item : *m_particles) {
			const QPointF position = sceneToViewport.map(QPointF(item.posX, item.posY));
			const QRectF particleRect(position.x() - radius, position.y() - radius, diameter, diameter);
			if (item.id < 8000) {
				const int colorGroup = item.id < 1500 ? 0 : item.id < 5000 ? 1 : item.id < 7500 ? 2 : 3;
				m_fixedColorPaths[colorGroup].addEllipse(particleRect);
				continue;
			}

			const QRgb color = QColor::fromRgbF(item.colorR, item.colorG, item.colorB).rgba();
			auto colorIndex = m_dynamicColorIndices.find(color);
			qsizetype index;
			if (colorIndex == m_dynamicColorIndices.end() || colorIndex.value().generation != m_colorGeneration) {
				index = static_cast<qsizetype>(m_activeDynamicColorPaths++);
				if (static_cast<std::size_t>(index) == m_dynamicColorPaths.size()) {
					m_dynamicColorPaths.emplace_back(color, QPainterPath{});
					m_dynamicColorPaths.back().second.setFillRule(Qt::WindingFill);
				} else {
					m_dynamicColorPaths[static_cast<std::size_t>(index)].first = color;
				}
				m_dynamicColorIndices.insert(color, ColorPathIndex { index, m_colorGeneration });
			} else {
				index = colorIndex.value().index;
			}
			m_dynamicColorPaths[static_cast<std::size_t>(index)].second.addEllipse(particleRect);
		}

		for (std::size_t group = 0; group < m_fixedColorPaths.size(); ++group) {
			painter->setBrush(fixedColors[group]);
			painter->drawPath(m_fixedColorPaths[group]);
		}

		// The dynamic-color batches are kept in first-seen order and drawn after
		// the contiguous fixed-color ID ranges, as in the previous renderer.
		for (std::size_t index = 0; index < m_activeDynamicColorPaths; ++index) {
			const auto& [color, path] = m_dynamicColorPaths[index];
			painter->setBrush(QColor::fromRgba(color));
			painter->drawPath(path);
		}
		painter->restore();

		if constexpr (performanceLoggingEnabled) {
			const qint64 paintNs = paintTimer.nsecsElapsed();
			if (paintNs >= 16000000) {
				logPerformanceEvent(QStringLiteral("slow_paint"),
				    QStringLiteral("paint_ms=%1 particles=%2 dynamic_color_batches=%3")
				        .arg(paintNs / 1000000.0, 0, 'f', 2)
				        .arg(m_particles->size())
				        .arg(m_activeDynamicColorPaths));
			}
		}
	}

private:
	struct ColorPathIndex {
		qsizetype index;
		quint32 generation;
	};

	const std::vector<Particle>* m_particles;
	std::array<QPainterPath, 4> m_fixedColorPaths;
	QHash<QRgb, ColorPathIndex> m_dynamicColorIndices;
	std::vector<std::pair<QRgb, QPainterPath>> m_dynamicColorPaths;
	std::size_t m_activeDynamicColorPaths = 0;
	quint32 m_colorGeneration = 0;
};

fluid::fluid(QWidget* parent) :
    ui(new Ui::fluid), m_paused(true), m_impl(std::make_unique<Impl>()), m_pointsPixmap(constants::maxwidth, constants::maxheight), m_pixmapSceneItem(nullptr)
{
	ui->setupUi(this);
	resetPerformanceLog();

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


	// Use a repeating timer so frame work does not add another 16 ms of idle
	// time after every simulation step.
	m_drawTimer.setTimerType(Qt::PreciseTimer);
	m_drawTimer.setInterval(16);
	connect(&m_drawTimer, &QTimer::timeout, this, &fluid::update);
	m_drawTimer.start();
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
		auto result = m_impl->f.simulate(constants::dt, performanceLoggingEnabled);
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
	float frameIntervalMs = 0.0f;
	QElapsedTimer frameWorkTimer;
	if constexpr (performanceLoggingEnabled) {
		frameIntervalMs = m_frameTimer.isValid() ? m_frameTimer.nsecsElapsed() * 0.000001f : 0.0f;
		m_frameTimer.restart();
		frameWorkTimer.start();
	}


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
	if constexpr (performanceLoggingEnabled) {
		const float frameWorkMs = frameWorkTimer.nsecsElapsed() * 0.000001f;
		if (frameIntervalMs >= 40.0f || frameWorkMs >= 40.0f) {
			logPerformanceEvent(QStringLiteral("slow_frame"),
			    QStringLiteral("frame=%1 interval_ms=%2 work_ms=%3 simulate_ms=%4 draw_request_ms=%5")
			        .arg(fps)
			        .arg(frameIntervalMs, 0, 'f', 2)
			        .arg(frameWorkMs, 0, 'f', 2)
			        .arg(simElapsed, 0, 'f', 2)
			        .arg(simDrawElapsed, 0, 'f', 2));
	}
	}

	m_timer.restart();

}


void fluid::drawImageToScene()
{
	if (m_pixmapSceneItem == nullptr) {
		m_pixmapSceneItem = m_pixmapScene->addPixmap(m_pointsPixmap);
		m_pixmapSceneItem->setCacheMode(QGraphicsItem::DeviceCoordinateCache);
		m_pixmapScene->setSceneRect(m_pointsPixmap.rect());
	}
}
