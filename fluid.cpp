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
	ui->m_graphView->scale(1, -1);	// invert Y axis
}

void fluid::setupScene()
{
	float res = 25;

	float tankHeight = 1.0 * constants::simheight;
	float tankWidth  = 1.0 * constants::simwidth;
	float h          = tankHeight / res;
	float density    = 1000.0f;

	float relWaterHeight = 0.8f;
	float relWaterWidth  = 0.6f;

	    // dam break

	    // compute number of particles

	float r = 0.5 * h; // particle radius w.r.t. cell size
	float dx    = 2.0 * r;
	float dy    = sqrt(3.0f) / 2.0f * dx;

	float numX         = floor((relWaterWidth * tankWidth - 2.0f * h - 2.0f * r) / dx);
	float numY         = floor((relWaterHeight * tankHeight - 2.0f * h - 2.0f * r) / dy);
	float maxParticles = numX * numY;

	// create fluid

	m_f = FlipFluid(density, tankWidth, tankHeight, h, r, maxParticles);

	// create particles

	m_f.numParticles = numX * numY;
	int p = 0;
	for (int i = 0; i < numX; i++) {
		for (int j = 0; j < numY; j++) {
			m_f.particlePos[p++] = h + r + dx * i + (j % 2 == 0 ? 0.0f : r);
			m_f.particlePos[p++] = h + r + dy * j;
		}
	}

	// setup grid cells for tank

	float n = m_f.fNumY;

	for (int i = 0; i < m_f.fNumX; i++) {
		for (int j = 0; j < m_f.fNumY; j++) {
			float s = 1.0f; // fluid
			if (i == 0 || i == m_f.fNumX - 1 || j == 0)
				s = 0.0f; // solid
			m_f.s[i * n + j] = s;
		}
	}


	//setup grid
	int nbGridItems = m_f.fNumX * m_f.fNumY;
	m_gridItems.reserve(nbGridItems);
	 if (constants::showGrid) {
		auto gridPen   = QPen(QColor(127, 127, 127));
		for (auto i = 0; i < m_f.fNumX; i++) {
			for (auto j = 0; j < m_f.fNumY; j++) {
				m_gridItems.push_back(m_scene->addRect(i * h,
				    j * h, h, h,
				    gridPen,
				    QBrush(QColor::fromRgbF(m_f.cellColor[3 * i], m_f.cellColor[3 * i + 1], m_f.cellColor[3 * i + 2]))));
			}
		}
	}

	//setup particles of water
	m_particleItems.reserve(m_f.numParticles);
	auto r2 = r * 2;
	 if (constants::showParticles) {
		for (auto i = 0; i < m_f.numParticles; i++) {
			auto particleColor = QColor::fromRgbF(m_f.particleColor[3 * i], m_f.particleColor[3 * i + 1], m_f.particleColor[3 * i + 2]);
			m_particleItems.push_back(m_scene->addEllipse(m_f.particlePos[2 * i],
			    m_f.particlePos[2 * i + 1],
			    r2, r2, QPen(particleColor), QBrush(particleColor)));
		}
	}

	setupObstacle(3.0, 2.0, true);
}

void fluid::setupObstacle(int x, int y, bool reset)
{
	float vx = 0.0;
	float vy = 0.0;

	if (!reset) {
		vx = (x - m_obstacleX) / m_scenedt;
		vy = (y - m_obstacleY) / m_scenedt;
	}

	m_obstacleX = x;
	m_obstacleY = y;
	auto r           = constants::obstacleRadius;
	auto n           = m_f.fNumY;
	auto cd          = sqrt(2) * m_f.h;

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

void fluid::simulate() {
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

void fluid::draw() {
	//grid
	auto p = 0;


	
	 if (constants::showGrid) {
		for (auto i = 0; i < m_f.fNumCells; i++) {
				m_gridItems[i]->setBrush(QBrush(QColor::fromRgbF(m_f.cellColor[3 * i], m_f.cellColor[3 * i + 1], m_f.cellColor[3 * i + 2])));
			}
		}
	

	//particles

	 //setup particles of water
	 m_particleItems.reserve(m_f.numParticles);
	 if (constants::showParticles) {
		 for (auto i = 0; i < m_f.numParticles; i++) {
			 auto particleColor = QColor::fromRgbF(m_f.particleColor[3 * i], m_f.particleColor[3 * i + 1], m_f.particleColor[3 * i + 2]);
			 m_particleItems.at(i)->setPos(m_f.particlePos[2 * i], m_f.particlePos[2 * i + 1]);
			 m_particleItems.at(i)->setBrush(QBrush(particleColor));
			 m_particleItems.at(i)->setPen(QPen(particleColor));
		 }
	 }

	 // obstacle

}

void fluid::requestAnimationFrame() {
}


void fluid::update() {
	QElapsedTimer timer;
	timer.start();

	simulate();
	draw();

	
	ui->label->setText(QString::number(timer.elapsed()) + " ms");
	//requestAnimationFrame();
	//callback to update
	//callback in 0ms
	QTimer::singleShot(0, this, &fluid::update);
}
