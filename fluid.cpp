#include "fluid.h"

#include <QtConcurrent>
#include <algorithm>
#include <execution>

#include <QElapsedTimer>

#include <hitPosition.h>

#include <calc.h>

#include "ui_fluid.h"


fluid::fluid(QWidget* parent) : ui(new Ui::fluid), m_frameNr(0), m_paused(true)
{
	ui->setupUi(this);
	m_scene = new QGraphicsScene(this);

	ui->m_graphView->setScene(m_scene);
	ui->m_graphView->scale(1.0, -1.0); // invert Y axis
}

void fluid::setupScene()
{
	m_scene->clear();
	m_particleItems.clear();
	m_gridItems.clear();

	double res = 20;

	double tankHeight = constants::simheight * constants::scale;
	double tankWidth  = constants::simwidth * constants::scale;
	double h          = tankHeight / res;
	double density    = 1000.0;

	double relWaterHeight = 0.8;
	double relWaterWidth  = 0.6;

	// dam break

	// compute number of particles

	double r  = 0.3 * h; // particle radius w.r.t. cell size
	double dx = 2.0 * r;
	double dy = sqrt(3.0) / 2.0 * dx;

	double numX         = floor((relWaterWidth * tankWidth - 2.0 * h - 2.0 * r) / dx);
	double numY         = floor((relWaterHeight * tankHeight - 2.0 * h - 2.0 * r) / dy);
	double maxParticles = numX * numY;

	// create fluid
	qCritical() << density << tankWidth << tankHeight << h << r << maxParticles;
	qCritical() << fixed << qSetRealNumberPrecision(20) << numX << numY << dx << dy << h + r + dx << h + r + dy;

	m_f = FlipFluid(density, tankWidth, tankHeight, h, r, maxParticles);

	// create particles
	auto r2          = r * 2.0;
	m_f.numParticles = numX * numY;
	int p            = 0;
	for (int i = 0; i < numX; i++) {
		for (int j = 0; j < numY; j++) {
			m_f.particlePos[p++] = h + r + dx * i + (j % 2 == 0 ? 0.0 : r);
			m_f.particlePos[p++] = h + r + dy * j;
		}
	}

	qCritical() << fixed << qSetRealNumberPrecision(20) << "create" << m_f.particlePos[0] << m_f.particlePos[1] << m_f.particleVel[0]
	            << m_f.particleVel[1];

	// setup grid cells for tank

	double n = m_f.fNumY;

	for (int i = 0; i < m_f.fNumX; i++) {
		for (int j = 0; j < m_f.fNumY; j++) {
			double s = 1.0; // fluid
			if (i == 0 || i == m_f.fNumX - 1 || j == 0)
				s = 0.0; // solid
			m_f.s[i * n + j] = s;
		}
	}

	double sizeRect = h;

	//setup grid
	int nbGridItems = m_f.fNumX * m_f.fNumY;
	m_gridItems.reserve(nbGridItems);

	auto pointSizeGrid = m_f.h;
	if (constants::showGrid) {
		auto gridPen = QPen(QColor(127, 127, 127));
		for (auto i = 0; i < m_f.fNumX; i++) {
			for (auto j = 0; j < m_f.fNumY; j++) {
				auto item = m_scene->addRect(0,
				    0,
				    pointSizeGrid,
				    pointSizeGrid,
				    gridPen,
				    QBrush(QColor::fromRgbF(m_f.cellColor[3 * i + j], m_f.cellColor[3 * i + j + 1], m_f.cellColor[3 * i + j + 2])));
				m_gridItems.push_back(item);
				item->setPos(i * pointSizeGrid, j * pointSizeGrid);
			}
		}
	}

	//setup particles of water
	m_particleItems.reserve(m_f.numParticles);
	if (constants::showParticles) {
		auto pointSize = 2.0 * r;
		qCritical() << pointSize;
		for (auto i = 0; i < m_f.numParticles; i++) {
			auto particleColor = QColor::fromRgbF(m_f.particleColor[3 * i], m_f.particleColor[3 * i + 1], m_f.particleColor[3 * i + 2]);
			m_particleItems.push_back(m_scene->addEllipse(0.0, 0.0, pointSize, pointSize, QPen(particleColor), QBrush(particleColor)));
			m_particleItems.at(i)->setPos(m_f.particlePos[2 * i] - r, m_f.particlePos[2 * i + 1] - r);
		}
	}

	setupObstacle(0.0, 0.0, true);
	m_obstacleItem = m_scene->addEllipse(0, 0, constants::obstacleRadius * 2.0, constants::obstacleRadius * 2.0, QPen(Qt::red), QBrush(Qt::red));
	m_obstacleItem->setFlag(QGraphicsItem::ItemIsMovable);
}

void fluid::setupObstacle(double x, double y, bool reset)
{
	double vx = 0.0;
	double vy = 0.0;

	if (!reset) {
		vx = (x - m_obstacleX) / constants::dt;
		vy = (y - m_obstacleY) / constants::dt;
	}

	m_obstacleX = x;
	m_obstacleY = y;
	auto r      = constants::obstacleRadius;
	auto n      = m_f.fNumY;
	auto cd     = sqrt(2.0) * m_f.h;

	for (auto i = 1; i < m_f.fNumX - 2; i++) {
		for (auto j = 1; j < m_f.fNumY - 2; j++) {
			m_f.s[i * n + j] = 1.0;

			auto dx = (i + 0.5) * m_f.h - x;
			auto dy = (j + 0.5) * m_f.h - y;

			if (dx * dx + dy * dy < r * r) {
				m_f.s[i * n + j]       = 0.0;
				m_f.u[i * n + j]       = vx;
				m_f.u[(i + 1) * n + j] = vx;
				m_f.v[i * n + j]       = vy;
				m_f.v[i * n + j + 1]   = vy;
			}
		}
	}

	//scene.showObstacle = true;
	ggobstacleVelX = vx;
	ggobstacleVelY = vy;
}

void fluid::simulate()
{
	if (!m_paused)
		m_f.simulate(constants::dt,
		    constants::gravity,
		    constants::flipRatio,
		    constants::numPressureIters,
		    constants::numParticleIters,
		    constants::overRelaxation,
		    constants::compensateDrift,
		    constants::separateParticles,
		    m_obstacleX,
		    m_obstacleY,
		    constants::obstacleRadius);

	m_frameNr++;
}

void fluid::draw()
{
	//grid
	auto p = 0;



	if (constants::showGrid) {
		for (auto i = 0; i < m_f.fNumCells; i++) {
			m_gridItems[i]->setPen(QPen(QColor(0, 0, 0)));
			m_gridItems[i]->setBrush(QBrush(QColor::fromRgbF(m_f.cellColor[3 * i], m_f.cellColor[3 * i + 1], m_f.cellColor[3 * i + 2])));
		}
	}


	//particles

	//setup particles of water
	if (constants::showParticles) {
		for (auto i = 0; i < m_f.numParticles; i++) {
			auto particleColor = QColor::fromRgbF(m_f.particleColor[3 * i], m_f.particleColor[3 * i + 1], m_f.particleColor[3 * i + 2]);

			m_particleItems.at(i)->setPos(m_f.particlePos[2 * i] - m_f.particleRadius, m_f.particlePos[2 * i + 1] - m_f.particleRadius);
			m_particleItems.at(i)->setBrush(QBrush(particleColor));
			m_particleItems.at(i)->setPen(QPen(particleColor));
		}
	}

	// obstacle
}

void fluid::requestAnimationFrame() {}


void fluid::update()
{
	QElapsedTimer timer;
	timer.start();

	setupObstacle(((double)m_obstacleItem->x() + (constants::obstacleRadius)), ((double)m_obstacleItem->y() + (constants::obstacleRadius)), false);

	//while (timer.elapsed() < 16.666f)
	simulate();
	draw();


	ui->label->setText(QString::number(timer.elapsed()) + " ms");
	ui->dsBox1->setValue(m_f.particlePos[0]);
	ui->dsBox2->setValue(m_f.particlePos[1]);
	ui->dsBox3->setValue(m_f.particleVel[0]);
	ui->dsBox4->setValue(m_f.particleVel[1]);
	ui->dsBox5->setValue(constants::scale);
	//ui->dsBox6->setValue(res);
	ui->dsBox7->setValue(m_f.h);
	ui->dsBox8->setValue(m_obstacleX);
	ui->dsBox9->setValue(m_obstacleY);
	ui->dsBox10->setValue(m_f.fNumY);
	ui->dsBox11->setValue(m_f.pNumY);

	auto h = 1.0 / m_f.fInvSpacing;
	auto r = m_f.particleRadius;

	auto minX = h + r;
	auto maxX = (m_f.fNumX - 1) * h - r;
	auto minY = h + r;
	auto maxY = (m_f.fNumY - 1) * h - r;

	ui->dsBox12->setValue(maxY);

	//requestAnimationFrame();
	//callback to update
	//callback in 0ms
	QTimer::singleShot(16, this, &fluid::update);
}
