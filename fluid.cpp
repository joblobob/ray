#include "fluid.h"

#include <QTimer>
#include <execution>

#include "ui_fluid.h"


fluid::fluid(QWidget* parent) : ui(new Ui::fluid), m_paused(true)
{
	ui->setupUi(this);
	m_scene = new QGraphicsScene(this);

	m_scene->setSceneRect(-50, -50, 550, 550);
	ui->m_graphView->setScene(m_scene);
	ui->m_graphView->scale(1.0, -1.0); // invert Y axis
}

void fluid::setupScene()
{
	m_scene->clear();
	m_particleItems.clear();
	m_gridItems.clear();

	double res = 50;

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
		auto gridPen         = QPen(QColor(127, 127, 127));
		int yesMyIndexAtHand = 0;
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
				item->setData(0, yesMyIndexAtHand);
				yesMyIndexAtHand++;
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
			m_particleItems.at(i)->setData(0, i);
		}
	}

	m_f.setupObstacle(0.0, 0.0, true);
	m_obstacleItem = m_scene->addEllipse(0, 0, constants::obstacleRadius * 2.0, constants::obstacleRadius * 2.0, QPen(Qt::red), QBrush(Qt::red));
	m_obstacleItem->setFlag(QGraphicsItem::ItemIsMovable);
}

void fluid::simulate()
{
	if (!m_paused) {
		auto result = m_f.simulate(constants::dt,
		    constants::gravity,
		    constants::flipRatio,
		    constants::numPressureIters,
		    constants::numParticleIters,
		    constants::overRelaxation,
		    constants::compensateDrift,
		    constants::separateParticles,
		    constants::obstacleRadius,
		    true);
		if (result.size() > 0) {
			ui->label_2->setText(result[0].message + ": " + QString::number(result[0].nsElapsed * 0.001) + " us");
			ui->label_3->setText(result[1].message + ": " + QString::number(result[1].nsElapsed * 0.001) + " us");
			ui->label_4->setText(result[2].message + ": " + QString::number(result[2].nsElapsed * 0.001) + " us");
			ui->label_5->setText(result[3].message + ": " + QString::number(result[3].nsElapsed * 0.001) + " us");
			ui->label_6->setText(result[4].message + ": " + QString::number(result[4].nsElapsed * 0.001) + " us");
			ui->label_7->setText(result[5].message + ": " + QString::number(result[5].nsElapsed * 0.001) + " us");
			ui->label_8->setText(result[6].message + ": " + QString::number(result[6].nsElapsed * 0.001) + " us");
			ui->label_9->setText(result[7].message + ": " + QString::number(result[7].nsElapsed * 0.001) + " us");
			ui->label_10->setText(result[8].message + ": " + QString::number(result[8].nsElapsed * 0.001) + " us");
		}
	}
}

void fluid::draw()
{
	//grid
	/* if (constants::showGrid) {
		auto setCellColor = [&](QGraphicsRectItem* item) {
			auto i = item->data(0).toInt();
			item->setBrush(QBrush(QColor::fromRgbF(m_f.cellColor[3 * i], m_f.cellColor[3 * i + 1], m_f.cellColor[3 * i + 2])));
		};

		std::for_each(std::execution::par_unseq, m_gridItems.begin(), m_gridItems.end(), setCellColor);
	}*/


	//particles
	//setup particles of water
	if (constants::showParticles) {
		/*auto setParticleColor = [&](QGraphicsEllipseItem* item) {
			auto i             = item->data(0).toInt();
			auto particleColor = QColor::fromRgbF(m_f.particleColor[3 * i], m_f.particleColor[3 * i + 1], m_f.particleColor[3 * i + 2]);
			item->setPos(m_f.particlePos[2 * i] - m_f.particleRadius, m_f.particlePos[2 * i + 1] - m_f.particleRadius);
			item->setBrush(QBrush(particleColor));
			item->setPen(QPen(particleColor));
		};

		std::for_each(std::execution::seq, m_particleItems.begin(), m_particleItems.end(), setParticleColor);
		*/

		for (auto i = 0; i < m_f.numParticles; i++) {
			auto particleColor = QColor::fromRgbF(m_f.particleColor[3 * i], m_f.particleColor[3 * i + 1], m_f.particleColor[3 * i + 2]);

			m_particleItems.at(i)->setPos(m_f.particlePos[2 * i] - m_f.particleRadius, m_f.particlePos[2 * i + 1] - m_f.particleRadius);
			m_particleItems.at(i)->setBrush(QBrush(particleColor));
			m_particleItems.at(i)->setPen(QPen(particleColor));
		}
	}
}



void fluid::update()
{
	m_timer.restart();

	m_f.setupObstacle(((double)m_obstacleItem->x() + (constants::obstacleRadius)), ((double)m_obstacleItem->y() + (constants::obstacleRadius)), false);

	//while (timer.elapsed() < 16.666f)
	simulate();
	auto elapsed = m_timer.nsecsElapsed() * 0.001;
	m_timer.restart();
	draw();

	ui->label->setText("Simulate: " + QString::number(elapsed) + " us" + " Draw: " + QString::number(m_timer.nsecsElapsed() * 0.001) + " us");


	//requestAnimationFrame();
	//callback to update
	//callback in 0ms
	QTimer::singleShot(16, this, &fluid::update);
}
