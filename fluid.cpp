#include "fluid.h"

#include <QGraphicsView>
#ifndef __EMSCRIPTEN__
#include <QOpenGLWidget>
#endif
#include <QCoreApplication>
#include <QDateTime>
#include <QFile>
#include <QPainter>
#include <QTextStream>
#include <QTimer>

#include <algorithm>
#include <array>
#include <cmath>
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

struct FixedParticleType {
	int count;
	QColor color;
};

constexpr std::size_t fixedParticleTypeCount = 4;
using FixedParticleTypes = std::array<FixedParticleType, fixedParticleTypeCount>;

const FixedParticleTypes& fixedParticleTypes()
{
	static const FixedParticleTypes types {
	    FixedParticleType { constants::rockParticleCount, QColor(Qt::darkGray) },
	    FixedParticleType { constants::coarseSandParticleCount, QColor(161, 103, 23) },
	    FixedParticleType { constants::fineSandParticleCount, QColor(191, 133, 53) },
	    FixedParticleType { constants::grassParticleCount, QColor(Qt::darkGreen) },
	};
	return types;
}

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
	float lastParticlePaintMs = 0.0f;
};

class ParticleOverlayItem final : public QGraphicsItem
{
public:
	explicit ParticleOverlayItem(const ParticleStore* particles, float* lastPaintMs) : m_particles(particles), m_lastPaintMs(lastPaintMs)
	{
		setZValue(1.0);
	}

	QRectF boundingRect() const override
	{
		return { 0.0, 0.0, constants::maxwidth, constants::maxheight };
	}

	void paint(QPainter* painter, const QStyleOptionGraphicsItem*, QWidget* widget) override
	{
		if (!constants::showParticles || m_particles == nullptr)
			return;
		QElapsedTimer paintTimer;
		const bool measurePaint = performanceLoggingEnabled || (++m_paintMeasureCounter % 4 == 0);
		if (measurePaint)
			paintTimer.start();

		// Rasterize into a reusable transparent image. Building a large path of
		// ellipses every frame is particularly expensive in the Wasm raster backend.
		const QTransform sceneToViewport = painter->worldTransform();
		const QPointF sceneOrigin = sceneToViewport.map(QPointF(0.0, 0.0));
		const QPointF screenSpacing = sceneToViewport.map(QPointF(constants::particleDiameter, 0.0)) - sceneOrigin;
		const qreal radius = std::hypot(screenSpacing.x(), screenSpacing.y()) * 0.5;
		const qreal devicePixelRatio = widget != nullptr ? widget->devicePixelRatioF() : 1.0;
		const QSize logicalSize = widget != nullptr ? widget->size() : QSize(painter->device()->width(), painter->device()->height());
		const QSize pixelSize(std::max(1, static_cast<int>(std::ceil(logicalSize.width() * devicePixelRatio))),
		    std::max(1, static_cast<int>(std::ceil(logicalSize.height() * devicePixelRatio))));
		if (m_particleImage.size() != pixelSize)
			m_particleImage = QImage(pixelSize, QImage::Format_ARGB32_Premultiplied);
		m_particleImage.setDevicePixelRatio(devicePixelRatio);
		m_particleImage.fill(Qt::transparent);
		const qreal pixelRadius = radius * devicePixelRatio;
		const qreal radiusSquared = pixelRadius * pixelRadius;
		const auto& particleTypes = fixedParticleTypes();

		for (std::size_t particleIndex = 0; particleIndex < m_particles->size(); ++particleIndex) {
			const int particleId = static_cast<int>(particleIndex);
			const QPointF position = sceneToViewport.map(QPointF(m_particles->posX[particleIndex], m_particles->posY[particleIndex])) * devicePixelRatio;
			int remainingTypeId = particleId;
			std::size_t colorGroup = 0;
			while (colorGroup < particleTypes.size() && remainingTypeId >= particleTypes[colorGroup].count) {
				remainingTypeId -= particleTypes[colorGroup].count;
				++colorGroup;
			}
			const QRgb color = colorGroup < particleTypes.size()
			    ? particleTypes[colorGroup].color.rgba()
			    : qRgba(qRound(std::clamp(m_particles->colorR[particleIndex], 0.0f, 1.0f) * 255.0f),
				  qRound(std::clamp(m_particles->colorG[particleIndex], 0.0f, 1.0f) * 255.0f),
				  qRound(std::clamp(m_particles->colorB[particleIndex], 0.0f, 1.0f) * 255.0f),
				  255);

			const int top = std::max(0, static_cast<int>(std::ceil(position.y() - pixelRadius - 0.5)));
			const int bottom = std::min(pixelSize.height() - 1, static_cast<int>(std::floor(position.y() + pixelRadius - 0.5)));
			for (int y = top; y <= bottom; ++y) {
				const qreal dy = y + 0.5 - position.y();
				const qreal remainingRadiusSquared = radiusSquared - dy * dy;
				if (remainingRadiusSquared < 0.0)
					continue;
				const qreal halfWidth = std::sqrt(remainingRadiusSquared);
				const int left = std::max(0, static_cast<int>(std::ceil(position.x() - halfWidth - 0.5)));
				const int right = std::min(pixelSize.width() - 1, static_cast<int>(std::floor(position.x() + halfWidth - 0.5)));
				if (left <= right) {
					auto* row = reinterpret_cast<QRgb*>(m_particleImage.scanLine(y));
					std::fill(row + left, row + right + 1, color);
				}
			}
		}
		painter->save();
		painter->resetTransform();
		painter->drawImage(QPointF(0.0, 0.0), m_particleImage);
		painter->restore();
		qint64 paintNs = 0;
		if (measurePaint) {
			paintNs = paintTimer.nsecsElapsed();
		}
		if (measurePaint && m_lastPaintMs != nullptr)
			*m_lastPaintMs = paintNs * 0.000001f;

		if constexpr (performanceLoggingEnabled) {
			if (paintNs >= 16000000) {
				logPerformanceEvent(QStringLiteral("slow_paint"),
				    QStringLiteral("paint_ms=%1 particles=%2")
				        .arg(paintNs / 1000000.0, 0, 'f', 2)
				        .arg(m_particles->size()));
			}
		}
	}

private:
	const ParticleStore* m_particles;
	float* m_lastPaintMs;
	QImage m_particleImage;
	quint32 m_paintMeasureCounter = 0;
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

	// create particles
	ParticleStore& particles = m_impl->f.particleMap;
	for (std::size_t i = 0; i < constants::numX; i++) {
		for (std::size_t j = 0; j < constants::numY; j++) {
			const float posX = constants::particleStartX + constants::particleSpacingX * i + (j % 2 == 0 ? 0.0f : constants::particleStaggerX);
			const float posY = constants::particleStartY + constants::particleRowSpacing * j;
			const std::size_t particleIndex = i * constants::numY + j;
			particles.posX[particleIndex] = posX;
			particles.posY[particleIndex] = posY;
			particles.velX[particleIndex] = 0.0f;
			particles.velY[particleIndex] = 0.0f;
			particles.colorR[particleIndex] = 0.0f;
			particles.colorG[particleIndex] = 0.0f;
			particles.colorB[particleIndex] = 1.0f;
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
		m_particleOverlayItem = new ParticleOverlayItem(&m_impl->f.particleMap, &m_impl->lastParticlePaintMs);
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
#ifdef __EMSCRIPTEN__
	const bool refreshStatus = (fps % 8) == 0;
	#else
	const bool refreshStatus = true;
	#endif
	float totalElapsed = refreshStatus ? m_timer.nsecsElapsed() * 0.000001f : 0.0f;
	float frameIntervalMs = 0.0f;
	QElapsedTimer frameWorkTimer;
	if constexpr (performanceLoggingEnabled) {
		frameIntervalMs = m_frameTimer.isValid() ? m_frameTimer.nsecsElapsed() * 0.000001f : 0.0f;
		m_frameTimer.restart();
		frameWorkTimer.start();
	}


	if (refreshStatus)
		m_simtimer.restart();

	setupObstacle(m_impl->f.gridCells,
	    ((float)m_obstacleItem->x() + ObstacleConstants::radius),
	    ((float)m_obstacleItem->y() + ObstacleConstants::radius),
	    m_impl->f.obstacle,
	    false);

	simulate();
	float simElapsed = refreshStatus ? m_simtimer.nsecsElapsed() * 0.000001f : 0.0f;
	if (refreshStatus)
		m_simtimer.restart();
	//draw();
#ifdef __EMSCRIPTEN__
	// The Wasm build uses Qt's raster viewport. Keep advancing the simulation
	// every timer tick, but avoid rebuilding and painting every particle at 60 Hz.
	if ((fps & 1) == 0)
		drawimage();
#else
	drawimage();
#endif

	float simDrawElapsed = refreshStatus ? m_simtimer.nsecsElapsed() * 0.000001f : 0.0f;
	fps++;

	if (refreshStatus) {
		std::string fpsLabel = std::format("Frame: {} Simulate: {:.2f} ms SimDraw: {:.2f} ms Paint: {:.2f} ms roundtrip: {:.2f} ms total {:.2f}",
		    fps,
		    simElapsed,
		    simDrawElapsed,
		    m_impl->lastParticlePaintMs,
		    totalElapsed,
		    simElapsed + simDrawElapsed + totalElapsed);
		ui->label->setText(fpsLabel.data());
	}
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
