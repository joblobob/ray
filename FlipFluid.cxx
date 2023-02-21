module;

#include <QDebug>
#include <algorithm>
#include <execution>
#include <vector>

#include <QElapsedTimer.h>

module FlipFluid;

import ParticleIntegration;
import PushParticles;
import ParticleCollision;
import ParticleDensity;
import ParticleIncompressibility;


FlipFluid::FlipFluid(double density, double width, double height, double spacing, double particleRadius, int maxParticles) :
    density(density),
    fNumX(floor(width / spacing) + 1),
    fNumY(floor(height / spacing) + 1),
    fNumCells(fNumX * fNumY),
    // fluid
    particleRadius(particleRadius),
    // particles
    maxParticles(maxParticles),
    particleRestDensity(0.0),
    pInvSpacing(1.0 / (2.2 * particleRadius)),
    pNumX(floor(width * pInvSpacing) + 1),
    pNumY(floor(height * pInvSpacing) + 1),
    pNumCells(pNumX * pNumY),
    pBorder { .maxX = pNumX, .maxY = pNumY },
    fBorder { .maxX = fNumX, .maxY = fNumY },
    particleMap(maxParticles),
    gridCells(fNumCells)
{
	h           = std::max(width / (double)fNumX, height / (double)fNumY);
	fInvSpacing = 1.0 / h;
}


void FlipFluid::setupObstacle(double x, double y, bool reset)
{
	double vx = 0.0;
	double vy = 0.0;

	if (!reset) {
		vx = (x - obstacleX) / constants::dt;
		vy = (y - obstacleY) / constants::dt;
	}

	obstacleX = x;
	obstacleY = y;
	auto r    = constants::obstacleRadius;
	auto n    = fNumY;
	auto cd   = sqrt(2.0) * h;

	for (auto i = 1; i < fNumX - 2; i++) {
		for (auto j = 1; j < fNumY - 2; j++) {
			gridCells[i * n + j].s = 1.0;
			auto dx                = (i + 0.5) * h - x;
			auto dy                = (j + 0.5) * h - y;
			if (dx * dx + dy * dy < r * r) {
				gridCells[i * n + j].s       = 0.0;
				gridCells[i * n + j].u       = vx;
				gridCells[(i + 1) * n + j].u = vx;
				gridCells[i * n + j].v       = vy;
				gridCells[i * n + j + 1].v   = vy;
			}
		}
	}

	//scene.showObstacle = true;
	obstacleVelX = vx;
	obstacleVelY = vy;
}


void FlipFluid::transferVelocitiesToParticles()
{
	for_each(std::execution::par_unseq, particleMap.begin(), particleMap.end(), [&](Particle& particle) { parseVelocitiesParticles(particle); });
}

void FlipFluid::transferVelocitiesToGrid()
{
	auto prepareCellType = [&](Cell& cell) {
		cell.prevU    = cell.u;
		cell.prevV    = cell.v;
		cell.du       = 0.0;
		cell.dv       = 0.0;
		cell.u        = 0.0;
		cell.v        = 0.0;
		cell.cellType = isVeryCloseToZero(cell.s) ? constants::CellType::Solid : constants::CellType::Air;
	};
	for_each(std::execution::seq, gridCells.begin(), gridCells.end(), prepareCellType);

	auto setCellType = [&](const Particle& p) {
		int cellNr     = cellNumber(p.posX, p.posY, fBorder, fInvSpacing);
		auto& cellType = gridCells[cellNr].cellType;
		if (cellType == constants::CellType::Air)
			cellType = constants::CellType::Fluid;
	};
	for_each(std::execution::seq, particleMap.begin(), particleMap.end(), setCellType);


	for_each(std::execution::seq, particleMap.begin(), particleMap.end(), [&](Particle& particle) { parseVelocitiesV(particle); });

	for_each(std::execution::seq, gridCells.begin(), gridCells.end(), [&](Cell& cell) { restoreSolidCellsV(cell); });
}

void FlipFluid::updateParticleColors()
{
	auto h1 = fInvSpacing;

	auto parseParticleColors = [&](Particle& particle) {
		auto s = 0.01;

		particle.colorR = std::clamp(particle.colorR - s, 0.0, 1.0);
		particle.colorG = std::clamp(particle.colorG - s, 0.0, 1.0);
		particle.colorB = std::clamp(particle.colorB + s, 0.0, 1.0);

		int cellNr = cellNumber({ particle.posX, particle.posY }, fBorder, fInvSpacing);

		auto d0 = particleRestDensity;

		if (d0 > 0.0) {
			const auto relDensity = gridCells[cellNr].particleDensity / d0;
			if (relDensity < 0.7) {
				const auto s    = 0.8;
				particle.colorR = s;
				particle.colorG = s;
				particle.colorB = 1.0;
			}
		}
	};
	for_each(std::execution::par_unseq, particleMap.begin(), particleMap.end(), parseParticleColors);
}

void FlipFluid::setSciColor(Cell& cell, double val, double minVal, double maxVal)
{
	val      = std::min(std::max(val, minVal), maxVal - 0.0001);
	auto d   = maxVal - minVal;
	val      = isVeryCloseToZero(d) ? 0.5 : (val - minVal) / d;
	double m = 0.25;
	int num  = floor(val / m);
	double s = (val - num * m) / m;
	double r, g, b;

	switch (num) {
		case 0:
			r = 0.0;
			g = s;
			b = 1.0;
			break;
		case 1:
			r = 0.0;
			g = 1.0;
			b = 1.0 - s;
			break;
		case 2:
			r = s;
			g = 1.0;
			b = 0.0;
			break;
		case 3:
			r = 1.0;
			g = 1.0 - s;
			b = 0.0;
			break;
	}

	cell.colorR = r;
	cell.colorG = g;
	cell.colorB = b;
}

void FlipFluid::updateCellColors()
{
	for_each(std::execution::seq, gridCells.begin(), gridCells.end(), [](Cell& c) {
		c.colorR = 0.0;
		c.colorG = 0.0;
		c.colorB = 0.0;
	});

	auto calcCellColor = [&](Cell& cell) {
		if (cell.cellType == constants::CellType::Solid) {
			cell.colorR = 0.5;
			cell.colorG = 0.5;
			cell.colorB = 0.5;
		} else if (cell.cellType == constants::CellType::Fluid) {
			auto d = cell.particleDensity;
			if (particleRestDensity > 0.0)
				d /= particleRestDensity;
			setSciColor(cell, d, 0.0, 2.0);
		}
	};
	for_each(std::execution::par_unseq, gridCells.begin(), gridCells.end(), calcCellColor);
}

std::vector<ExecutionLog> FlipFluid::simulate(double dt,
    double flipRatio,
    int numPressureIters,
    int numParticleIters,
    double overRelaxation,
    double obstacleRadius,
    bool instrument)
{
	std::vector<ExecutionLog> log;



	QElapsedTimer timer;

	if (instrument)
		timer.start();

	integrateParticles(particleMap);
	if (instrument) {
		log.push_back(ExecutionLog { "integrateParticles", timer.nsecsElapsed() });
		timer.restart();
	}
	pushParticlesApart(particleMap, numParticleIters, pNumCells, pNumX, pNumY, maxParticles, pBorder, pInvSpacing, particleRadius);
	if (instrument) {
		log.push_back(ExecutionLog { "pushParticlesApart", timer.nsecsElapsed() });
		timer.restart();
	}
	handleParticleCollisions(particleMap, fNumX, fNumY, fInvSpacing, obstacleX, obstacleY, obstacleVelX, obstacleVelY, obstacleRadius, particleRadius);
	if (instrument) {
		log.push_back(ExecutionLog { "handleParticleCollisions", timer.nsecsElapsed() });
		timer.restart();
	}
	transferVelocitiesToGrid();
	if (instrument) {
		log.push_back(ExecutionLog { "transferVelocitiesTrue", timer.nsecsElapsed() });
		timer.restart();
	}
	updateParticleDensity(particleMap, gridCells, fNumX, fNumY, fInvSpacing, h);
	if (instrument) {
		log.push_back(ExecutionLog { "updateParticleDensity", timer.nsecsElapsed() });
		timer.restart();
	}
	if (isVeryCloseToZero(particleRestDensity)) {
		calculateParticleRestDensity(gridCells, particleRestDensity);
	}
	solveIncompressibility(gridCells, numPressureIters, fNumY, density, h, dt, overRelaxation, particleRestDensity);
	if (instrument) {
		log.push_back(ExecutionLog { "solveIncompressibility", timer.nsecsElapsed() });
		timer.restart();
	}
	transferVelocitiesToParticles();
	if (instrument) {
		log.push_back(ExecutionLog { "transferVelocitiesFalse", timer.nsecsElapsed() });
		timer.restart();
	}


	updateParticleColors();
	if (instrument) {
		log.push_back(ExecutionLog { "updateParticleColors", timer.nsecsElapsed() });
		timer.restart();
	}
	updateCellColors();
	if (instrument) {
		log.push_back(ExecutionLog { "updateCellColors", timer.nsecsElapsed() });
		timer.restart();
	}
	return log;
}

void FlipFluid::parseVelocitiesParticles(Particle& particle)
{
	double h2 = 0.5 * h;
	auto n    = fNumY;

	auto dx = 0.0;
	auto dy = h2;

	const auto x = std::clamp(particle.posX, h, (fNumX - 1) * h);
	const auto y = std::clamp(particle.posY, h, (fNumY - 1) * h);

	int x0  = std::min((int)floor((x - dx) * fInvSpacing), fNumX - 2);
	auto tx = ((x - dx) - x0 * h) * fInvSpacing;
	int x1  = std::min(x0 + 1, fNumX - 2);

	int y0  = std::min((int)floor((y - dy) * fInvSpacing), fNumY - 2);
	auto ty = ((y - dy) - y0 * h) * fInvSpacing;
	int y1  = std::min(y0 + 1, fNumY - 2);

	auto sx = 1.0 - tx;
	auto sy = 1.0 - ty;

	auto d0 = sx * sy;
	auto d1 = tx * sy;
	auto d2 = tx * ty;
	auto d3 = sx * ty;

	int nr0 = x0 * n + y0;
	int nr1 = x1 * n + y0;
	int nr2 = x1 * n + y1;
	int nr3 = x0 * n + y1;

	Cell& cell0 = gridCells[nr0];
	Cell& cell1 = gridCells[nr1];
	Cell& cell2 = gridCells[nr2];
	Cell& cell3 = gridCells[nr3];

	// specific

	int offset             = n;
	const Cell cellOffset0 = gridCells[nr0 - offset];
	const Cell cellOffset1 = gridCells[nr1 - offset];
	const Cell cellOffset2 = gridCells[nr2 - offset];
	const Cell cellOffset3 = gridCells[nr3 - offset];

	auto valid0 = cell0.cellType != constants::CellType::Air || cellOffset0.cellType != constants::CellType::Air ? 1.0 : 0.0;
	auto valid1 = cell1.cellType != constants::CellType::Air || cellOffset1.cellType != constants::CellType::Air ? 1.0 : 0.0;
	auto valid2 = cell2.cellType != constants::CellType::Air || cellOffset2.cellType != constants::CellType::Air ? 1.0 : 0.0;
	auto valid3 = cell3.cellType != constants::CellType::Air || cellOffset3.cellType != constants::CellType::Air ? 1.0 : 0.0;

	auto v = particle.velX;
	auto d = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

	if (d > 0.0) {
		auto picV = (valid0 * d0 * cell0.u + valid1 * d1 * cell1.u + valid2 * d2 * cell2.u + valid3 * d3 * cell3.u) / d;
		auto corr = (valid0 * d0 * (cell0.u - cell0.prevU) + valid1 * d1 * (cell1.u - cell1.prevU) + valid2 * d2 * (cell2.u - cell2.prevU) +
		                valid3 * d3 * (cell3.u - cell3.prevU)) /
		            d;
		auto flipV = v + corr;

		auto newVelValue = (1.0 - constants::flipRatio) * picV + constants::flipRatio * flipV;
		particle.velX    = newVelValue;
	}

	// 2eme passe
	dx = h2;
	dy = 0.0;

	x0 = std::min((int)floor((x - dx) * fInvSpacing), fNumX - 2);
	tx = ((x - dx) - x0 * h) * fInvSpacing;
	x1 = std::min(x0 + 1, fNumX - 2);

	y0 = std::min((int)floor((y - dy) * fInvSpacing), fNumY - 2);
	ty = ((y - dy) - y0 * h) * fInvSpacing;
	y1 = std::min(y0 + 1, fNumY - 2);

	sx = 1.0 - tx;
	sy = 1.0 - ty;

	d0 = sx * sy;
	d1 = tx * sy;
	d2 = tx * ty;
	d3 = sx * ty;

	nr0 = x0 * n + y0;
	nr1 = x1 * n + y0;
	nr2 = x1 * n + y1;
	nr3 = x0 * n + y1;

	Cell& cell0v = gridCells[nr0];
	Cell& cell1v = gridCells[nr1];
	Cell& cell2v = gridCells[nr2];
	Cell& cell3v = gridCells[nr3];

	// specific


	offset                  = 1;
	const Cell cellOffset0v = gridCells[nr0 - offset];
	const Cell cellOffset1v = gridCells[nr1 - offset];
	const Cell cellOffset2v = gridCells[nr2 - offset];
	const Cell cellOffset3v = gridCells[nr3 - offset];

	valid0 = cell0.cellType != constants::CellType::Air || cellOffset0v.cellType != constants::CellType::Air ? 1.0 : 0.0;
	valid1 = cell1.cellType != constants::CellType::Air || cellOffset1v.cellType != constants::CellType::Air ? 1.0 : 0.0;
	valid2 = cell2.cellType != constants::CellType::Air || cellOffset2v.cellType != constants::CellType::Air ? 1.0 : 0.0;
	valid3 = cell3.cellType != constants::CellType::Air || cellOffset3v.cellType != constants::CellType::Air ? 1.0 : 0.0;

	v = particle.velY;
	d = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

	if (d > 0.0) {
		auto picV = (valid0 * d0 * cell0v.v + valid1 * d1 * cell1v.v + valid2 * d2 * cell2v.v + valid3 * d3 * cell3v.v) / d;
		auto corr = (valid0 * d0 * (cell0v.v - cell0v.prevV) + valid1 * d1 * (cell1v.v - cell1v.prevV) + valid2 * d2 * (cell2v.v - cell2v.prevV) +
		                valid3 * d3 * (cell3v.v - cell3v.prevV)) /
		            d;
		auto flipV = v + corr;

		auto newVelValue = (1.0 - constants::flipRatio) * picV + constants::flipRatio * flipV;
		particle.velY    = newVelValue;
	}
}

void FlipFluid::parseVelocitiesV(Particle& particle)
{
	double h2 = 0.5 * h;
	auto n    = fNumY;

	auto dx = 0.0;
	auto dy = h2;

	const auto x = std::clamp(particle.posX, h, (fNumX - 1) * h);
	const auto y = std::clamp(particle.posY, h, (fNumY - 1) * h);

	int x0  = std::min((int)floor((x - dx) * fInvSpacing), fNumX - 2);
	auto tx = ((x - dx) - x0 * h) * fInvSpacing;
	int x1  = std::min(x0 + 1, fNumX - 2);

	int y0  = std::min((int)floor((y - dy) * fInvSpacing), fNumY - 2);
	auto ty = ((y - dy) - y0 * h) * fInvSpacing;
	int y1  = std::min(y0 + 1, fNumY - 2);

	auto sx = 1.0 - tx;
	auto sy = 1.0 - ty;

	auto d0 = sx * sy;
	auto d1 = tx * sy;
	auto d2 = tx * ty;
	auto d3 = sx * ty;

	int nr0 = x0 * n + y0;
	int nr1 = x1 * n + y0;
	int nr2 = x1 * n + y1;
	int nr3 = x0 * n + y1;

	Cell& cell0 = gridCells[nr0];
	Cell& cell1 = gridCells[nr1];
	Cell& cell2 = gridCells[nr2];
	Cell& cell3 = gridCells[nr3];

	//specific

	setVelComponent(cell0.u, cell0.du, particle.velX, d0);
	setVelComponent(cell1.u, cell1.du, particle.velX, d1);
	setVelComponent(cell2.u, cell2.du, particle.velX, d2);
	setVelComponent(cell3.u, cell3.du, particle.velX, d3);


	dx = h2;
	dy = 0.0;

	x0 = std::min((int)floor((x - dx) * fInvSpacing), fNumX - 2);
	tx = ((x - dx) - x0 * h) * fInvSpacing;
	x1 = std::min(x0 + 1, fNumX - 2);

	y0 = std::min((int)floor((y - dy) * fInvSpacing), fNumY - 2);
	ty = ((y - dy) - y0 * h) * fInvSpacing;
	y1 = std::min(y0 + 1, fNumY - 2);

	sx = 1.0 - tx;
	sy = 1.0 - ty;

	d0 = sx * sy;
	d1 = tx * sy;
	d2 = tx * ty;
	d3 = sx * ty;

	nr0 = x0 * n + y0;
	nr1 = x1 * n + y0;
	nr2 = x1 * n + y1;
	nr3 = x0 * n + y1;

	Cell& cell0v = gridCells[nr0];
	Cell& cell1v = gridCells[nr1];
	Cell& cell2v = gridCells[nr2];
	Cell& cell3v = gridCells[nr3];

	setVelComponent(cell0v.v, cell0v.dv, particle.velY, d0);
	setVelComponent(cell1v.v, cell1v.dv, particle.velY, d1);
	setVelComponent(cell2v.v, cell2v.dv, particle.velY, d2);
	setVelComponent(cell3v.v, cell3v.dv, particle.velY, d3);
}

void FlipFluid::setVelComponent(double& f, double& d, double pv, double delta)
{
	f += pv * delta;
	d += delta;
}


void FlipFluid::restoreSolidCellsV(Cell& cell)
{
	if (cell.du > 0.0)
		cell.u /= cell.du;

	if (cell.dv > 0.0)
		cell.v /= cell.dv;

	auto solid = cell.cellType == constants::CellType::Solid;
	if (solid || (cell.cellNumX > 0 && gridCells[(cell.cellNumX - 1) * fNumY + cell.cellNumY].cellType == constants::CellType::Solid))
		cell.u = cell.prevU;
	if (solid || (cell.cellNumY > 0 && gridCells[cell.cellNumX * fNumY + cell.cellNumY - 1].cellType == constants::CellType::Solid))
		cell.v = cell.prevV;
}