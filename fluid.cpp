#include "fluid.h"

#include <QtConcurrent>
#include <algorithm>
#include <execution>

#include <QElapsedTimer>

#include <hitPosition.h>

#include <calc.h>

#include "ui_fluid.h"


fluid::fluid(QWidget* parent) : ui(new Ui::fluid), frameNr(0), paused(true)
{
	ui->setupUi(this);
	scene = new QGraphicsScene(this);

	ui->m_graphView->setScene(scene);

}

void fluid::setupScene()
{
	float res = 100;

	float tankHeight = 1.0 * constants::simheight;
	float tankWidth  = 1.0 * constants::simwidth;
	float h          = tankHeight / res;
	float density    = 1000.0;

	float relWaterHeight = 0.8;
	float relWaterWidth  = 0.6;

	    // dam break

	    // compute number of particles

	float r = 0.3 * h; // particle radius w.r.t. cell size
	float dx    = 2.0 * r;
	float dy    = sqrt(3.0) / 2.0 * dx;

	float numX         = floor((relWaterWidth * tankWidth - 2.0 * h - 2.0 * r) / dx);
	float numY         = floor((relWaterHeight * tankHeight - 2.0 * h - 2.0 * r) / dy);
	float maxParticles = numX * numY;

	// create fluid

	f = FlipFluid(density, tankWidth, tankHeight, h, r, maxParticles);

	// create particles

	f.numParticles = numX * numY;
	int p = 0;
	for (int i = 0; i < numX; i++) {
		for (int j = 0; j < numY; j++) {
			f.particlePos[p++] = h + r + dx * i + (j % 2 == 0 ? 0.0 : r);
			f.particlePos[p++] = h + r + dy * j;
		}
	}

	// setup grid cells for tank

	float n = f.fNumY;

	for (int i = 0; i < f.fNumX; i++) {
		for (int j = 0; j < f.fNumY; j++) {
			float s = 1.0; // fluid
			if (i == 0 || i == f.fNumX - 1 || j == 0)
				s = 0.0; // solid
			f.s[i * n + j] = s;
		}
	}

	//setup grid
	 if (constants::showGrid) {
	 auto pointSize = 0.9 * h / constants::simwidth * ui->m_graphView->width();
		for (auto i = 0; i < f.fNumX; i++) {
			for (auto j = 0; j < f.fNumY; j++) {
				scene->addRect(
				    (i + 0.5) * f.h, (j + 0.5) * f.h, 5, 5, QPen(QColor::fromRgbF(f.cellColor[3 * i], f.cellColor[3 * i + 1], f.cellColor[3 * i + 2])), Qt::SolidPattern);
			}
		}
	}

	//setup particles of water
	 if (constants::showParticles) {
		for (auto i = 0; i < f.numParticles; i++) {
			auto particlePen = QPen(QColor::fromRgbF(f.particleColor[3 * i], f.particleColor[3 * i + 1], f.particleColor[3 * i + 2]));
			scene->addEllipse(f.particlePos[2 * i], f.particlePos[2 * i + 1], r, r, particlePen, Qt::SolidPattern);
		}
	}

	setupObstacle(3.0, 2.0, true);
}

void fluid::setupObstacle(int x, int y, bool reset)
{
	float vx = 0.0;
	float vy = 0.0;

	if (!reset) {
		vx = (x - obstacleX) / scenedt;
		vy = (y - obstacleY) / scenedt;
	}

	obstacleX = x;
	obstacleY = y;
	auto r           = constants::obstacleRadius;
	auto n           = f.fNumY;
	auto cd          = sqrt(2) * f.h;

	for (auto i = 1; i < f.fNumX - 2; i++) {
		for (auto j = 1; j < f.fNumY - 2; j++) {
			f.s[i * n + j] = 1.0;

			auto dx = (i + 0.5) * f.h - x;
			auto dy = (j + 0.5) * f.h - y;

			if (dx * dx + dy * dy < r * r) {
				f.s[i * n + j]       = 0.0;
				f.u[i * n + j]       = vx;
				f.u[(i + 1) * n + j] = vx;
				f.v[i * n + j]       = vy;
				f.v[i * n + j + 1]   = vy;
			}
		}
	}

	//scene.showObstacle = true;
	ggobstacleVelX = vx;
	ggobstacleVelY = vy;
}

void fluid::simulate() {
	if (!paused)
		f.simulate(constants::dt,
		    constants::gravity,
		    constants::flipRatio,
		    constants::numPressureIters,
		    constants::numParticleIters,
		    constants::overRelaxation,
		    constants::compensateDrift,
		    constants::separateParticles,
		    obstacleX,
		    obstacleY,
		    constants::obstacleRadius);

	frameNr++;
}

void fluid::draw() {
	//grid
	auto p = 0;


	
	/* if (constants::showGrid) {
		auto pointSize = 0.9 * f.h / constants::simwidth * ui->m_graphView->width();
		QPen outlinePen(Qt::black);
		for (auto i = 0; i < f.fNumX; i++) {
			for (auto j = 0; j < f.fNumY; j++) {
				auto item = ui->m_graphView->itemAt((i + 0.5) * f.h, (j + 0.5) * f.h);
				if (item != nullptr) {
					auto rect = qgraphicsitem_cast<QGraphicsRectItem*>(item);
					if (rect != nullptr)
						rect->setPen(QPen(QColor::fromRgbF(f.cellColor[3 * i], f.cellColor[3 * i + 1], f.cellColor[3 * i + 2])));
				}
			}
		}
	}*/

	//particles

	auto pointSize = 2.0 * f.particleRadius / constants::simwidth * ui->m_graphView->width();

}

void fluid::requestAnimationFrame() {
}


void fluid::update() {
	QElapsedTimer timer;
	timer.start();

	simulate();
	draw();

	
	ui->label->setText(QString::number(timer.elapsed()) + " ms " + QString::number(f.particlePos[50]));
	//requestAnimationFrame();
	//callback to update
	//callback in 0ms
	QTimer::singleShot(0, this, &fluid::update);
}
