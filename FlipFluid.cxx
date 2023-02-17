module;

#include <QDebug>
#include <algorithm>
#include <execution>
#include <vector>

#include <QElapsedTimer.h>

module FlipFluid;


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
    cellParticleIds(maxParticles),
    pBorder { .maxX = pNumX, .maxY = pNumY },
    fBorder { .maxX = fNumX, .maxY = fNumY },
    particleMap(maxParticles),
    particleCells(pNumCells + 1),
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

void FlipFluid::integrateParticles()
{
	auto integration = [](Particle& particle) {
		particle.velY += constants::integ; //dt * gravity
		particle.posX += particle.velX * constants::dt;
		particle.posY += particle.velY * constants::dt;
	};
	std::for_each(std::execution::seq, particleMap.begin(), particleMap.end(), integration);
}

void FlipFluid::pushParticlesApart(int numIters)
{
	for_each(std::execution::seq, particleCells.begin(), particleCells.end(), [](ParticleInCells& p) { p.numCellParticles = 0; });

	std::vector<int> cellNumbers(maxParticles);

	// count particles per cell
	auto countParticlesInCell = [&](const Particle& particle) {
		const int cellnum        = cellNumber(particle.posX, particle.posY, pBorder, pInvSpacing);
		cellNumbers[particle.id] = cellnum;
		particleCells[cellnum].numCellParticles++;
	};
	for_each(std::execution::seq, particleMap.begin(), particleMap.end(), countParticlesInCell);

	// partial sums
	auto first = 0;

	//count first cell
	auto countFirstCell = [&](ParticleInCells& pIncell) {
		first += pIncell.numCellParticles;
		pIncell.firstCellParticle = first;
	};
	for_each(std::execution::seq, particleCells.begin(), particleCells.end(), countFirstCell);


	particleCells[pNumCells].firstCellParticle = first; // guard

	// fill particles into cells
	auto fillParticleIntoCell = [&](const Particle& p) {
		cellParticleIds[particleCells[cellNumbers[p.id]].firstCellParticle--] = p.id;
	};
	for_each(std::execution::seq, particleMap.begin(), particleMap.end(), fillParticleIntoCell);


	// push particles apart
	const auto minDist  = 2.0 * particleRadius;
	const auto minDist2 = minDist * minDist;

	auto pushParticles = [&](auto& particle) {
		const auto px = particle.posX;
		const auto py = particle.posY;

		const int pxi = floor(px * pInvSpacing);
		const int pyi = floor(py * pInvSpacing);
		const int x0  = std::max(pxi - 1, 0);
		const int y0  = std::max(pyi - 1, 0);
		const int x1  = std::min(pxi + 1, pNumX - 1);
		const int y1  = std::min(pyi + 1, pNumY - 1);

		for (int xi = x0; xi <= x1; xi++) {
			for (int yi = y0; yi <= y1; yi++) {
				const int cellNr = xi * pNumY + yi;
				const int first  = particleCells[cellNr].firstCellParticle;
				const int last   = particleCells[cellNr + 1].firstCellParticle;
				for (int j = first; j < last; j++) {
					int id = cellParticleIds[j];
					if (id != particle.id) {
						auto& particleAtId { particleMap[id] };
						const auto qx = particleAtId.posX;
						const auto qy = particleAtId.posY;

						auto dx       = qx - px;
						auto dy       = qy - py;
						const auto d2 = dx * dx + dy * dy;
						if (!(d2 > minDist2 || isVeryCloseToZero(d2))) {
							const auto d = sqrt(d2);
							const auto s = 0.5 * (minDist - d) / d;
							dx *= s;
							dy *= s;
							particle.posX -= dx;
							particle.posY -= dy;
							particleAtId.posX += dx;
							particleAtId.posY += dy;

							// diffuse colors

							calcColor(particle.colorR, particleAtId.colorR);
							calcColor(particle.colorG, particleAtId.colorG);
							calcColor(particle.colorB, particleAtId.colorB);
						}
					}
				}
			}
		}
	};

	for (int iter = 0; iter < numIters; iter++) {
		for_each(std::execution::par_unseq, particleMap.begin(), particleMap.end(), pushParticles);
	}
}

constexpr void FlipFluid::calcColor(double& p1Color, double& p2Color)
{
	const auto color = (p1Color + p2Color) * 0.5;
	p1Color          = p1Color + (color - p1Color) * constants::colorDiffusionCoeff;
	p2Color          = p2Color + (color - p2Color) * constants::colorDiffusionCoeff;
}

void FlipFluid::handleParticleCollisions(double obstacleX, double obstacleY, double obstacleRadius)
{
	const auto h        = 1.0 / fInvSpacing;
	const auto r        = particleRadius;
	const auto minDist  = obstacleRadius + r;
	const auto minDist2 = minDist * minDist;

	const auto minX = h + r;
	const auto maxX = (fNumX - 1) * h - r;
	const auto minY = h + r;
	const auto maxY = (fNumY - 1) * h - r;

	auto particleCollision = [&](auto& particle) {
		auto x = particle.posX;
		auto y = particle.posY;

		const auto dx = x - obstacleX;
		const auto dy = y - obstacleY;
		const auto d2 = dx * dx + dy * dy;

		// obstacle collision

		if (d2 < minDist2) {
			particle.velX = obstacleVelX;
			particle.velY = obstacleVelY;
		}
		// wall collisions

		if (x < minX) {
			x             = minX;
			particle.velX = 0.0;
		} else if (x > maxX) {
			x             = maxX;
			particle.velX = 0.0;
		}
		if (y < minY) {
			y             = minY;
			particle.velY = 0.0;
		} else if (y > maxY) {
			y             = maxY;
			particle.velY = 0.0;
		}
		particle.posX = x;
		particle.posY = y;
	};
	for_each(std::execution::seq, particleMap.begin(), particleMap.end(), particleCollision);
}

void FlipFluid::updateParticleDensity()
{
	const int n     = fNumY;
	const auto h1   = fInvSpacing;
	const double h2 = 0.5 * h;

	for_each(std::execution::seq, gridCells.begin(), gridCells.end(), [](Cell& c) { c.particleDensity = 0.0; });


	auto solveParticleDensity = [&](const Particle& p) {
		const auto x = std::clamp(p.posX, h, (fNumX - 1) * h);
		const auto y = std::clamp(p.posY, h, (fNumY - 1) * h);

		const int x0  = floor((x - h2) * h1);
		const auto tx = ((x - h2) - x0 * h) * h1;
		const int x1  = std::min(x0 + 1, fNumX - 2);

		const int y0  = floor((y - h2) * h1);
		const auto ty = ((y - h2) - y0 * h) * h1;
		const int y1  = std::min(y0 + 1, fNumY - 2);

		const auto sx = 1.0 - tx;
		const auto sy = 1.0 - ty;

		if (x0 < fNumX && y0 < fNumY)
			gridCells[x0 * n + y0].particleDensity += sx * sy;
		if (x1 < fNumX && y0 < fNumY)
			gridCells[x1 * n + y0].particleDensity += tx * sy;
		if (x1 < fNumX && y1 < fNumY)
			gridCells[x1 * n + y1].particleDensity += tx * ty;
		if (x0 < fNumX && y1 < fNumY)
			gridCells[x0 * n + y1].particleDensity += sx * ty;
	};
	for_each(std::execution::par_unseq, particleMap.begin(), particleMap.end(), solveParticleDensity);


	if (isVeryCloseToZero(particleRestDensity)) {
		double sum        = 0.0;
		int numFluidCells = 0;

		auto calculateSum = [&](Cell& cell) {
			if (cell.cellType == constants::CellType::Fluid) {
				sum += cell.particleDensity;
				numFluidCells++;
			}
		};
		for_each(std::execution::seq, gridCells.begin(), gridCells.end(), calculateSum);

		if (numFluidCells > 0)
			particleRestDensity = sum / numFluidCells;
	}
}

void FlipFluid::transferVelocitiesToParticles()
{
	auto n    = fNumY;
	auto h1   = fInvSpacing;
	double h2 = 0.5 * h;


	for (int component = 0; component < 2; component++) {
		auto dx = component == 0 ? 0.0 : h2;
		auto dy = component == 0 ? h2 : 0.0;

		auto parseVelocities = [&](Particle& particle) {
			const auto x = std::clamp(particle.posX, h, (fNumX - 1) * h);
			const auto y = std::clamp(particle.posY, h, (fNumY - 1) * h);

			int x0  = std::min((int)floor((x - dx) * h1), fNumX - 2);
			auto tx = ((x - dx) - x0 * h) * h1;
			int x1  = std::min(x0 + 1, fNumX - 2);

			int y0  = std::min((int)floor((y - dy) * h1), fNumY - 2);
			auto ty = ((y - dy) - y0 * h) * h1;
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

			double& cf0     = component == 0 ? cell0.u : cell0.v;
			double& cprevF0 = component == 0 ? cell0.prevU : cell0.prevV;

			double& cf1     = component == 0 ? cell1.u : cell1.v;
			double& cprevF1 = component == 0 ? cell1.prevU : cell1.prevV;

			double& cf2     = component == 0 ? cell2.u : cell2.v;
			double& cprevF2 = component == 0 ? cell2.prevU : cell2.prevV;

			double& cf3     = component == 0 ? cell3.u : cell3.v;
			double& cprevF3 = component == 0 ? cell3.prevU : cell3.prevV;



			int offset             = component == 0 ? n : 1;
			const Cell cellOffset0 = gridCells[nr0 - offset];
			const Cell cellOffset1 = gridCells[nr1 - offset];
			const Cell cellOffset2 = gridCells[nr2 - offset];
			const Cell cellOffset3 = gridCells[nr3 - offset];

			auto valid0 = cell0.cellType != constants::CellType::Air || cellOffset0.cellType != constants::CellType::Air ? 1.0 : 0.0;
			auto valid1 = cell1.cellType != constants::CellType::Air || cellOffset1.cellType != constants::CellType::Air ? 1.0 : 0.0;
			auto valid2 = cell2.cellType != constants::CellType::Air || cellOffset2.cellType != constants::CellType::Air ? 1.0 : 0.0;
			auto valid3 = cell3.cellType != constants::CellType::Air || cellOffset3.cellType != constants::CellType::Air ? 1.0 : 0.0;

			auto v = component == 0 ? particle.velX : particle.velY;
			auto d = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

			if (d > 0.0) {
				auto picV  = (valid0 * d0 * cf0 + valid1 * d1 * cf1 + valid2 * d2 * cf2 + valid3 * d3 * cf3) / d;
				auto corr  = (valid0 * d0 * (cf0 - cprevF0) + valid1 * d1 * (cf1 - cprevF1) + valid2 * d2 * (cf2 - cprevF2) + valid3 * d3 * (cf3 - cprevF3)) / d;
				auto flipV = v + corr;

				auto newVelValue = (1.0 - constants::flipRatio) * picV + constants::flipRatio * flipV;
				if (component == 0)
					particle.velX = newVelValue;
				else
					particle.velY = newVelValue;
			}
		};
		for_each(std::execution::par_unseq, particleMap.begin(), particleMap.end(), parseVelocities);
	}
}

void FlipFluid::transferVelocitiesToGrid()
{
	auto n    = fNumY;
	auto h1   = fInvSpacing;
	double h2 = 0.5 * h;

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


	auto parseVelocities = [&](Particle& particle) {
		parseVelocitiesV(particle);
	};
	for_each(std::execution::seq, particleMap.begin(), particleMap.end(), parseVelocities);


	auto restoreSolidCells = [&](Cell& cell) {
		restoreSolidCellsV(cell);
	};
	for_each(std::execution::seq, gridCells.begin(), gridCells.end(), restoreSolidCells);
}


void FlipFluid::solveIncompressibility(int numIters, double dt, double overRelaxation, bool compensateDrift)
{
	auto n  = fNumY;
	auto cp = density * h / dt;

	auto resetPrev = [&](Cell& cell) {
		cell.prevU = cell.u;
		cell.prevV = cell.v;
	};
	for_each(std::execution::seq, gridCells.begin(), gridCells.end(), resetPrev);

	for (int iter = 0; iter < numIters; iter++) {
		auto parseIncompressibility = [&](Cell& cell) {
			if (cell.cellType == constants::CellType::Fluid) {
				int i      = cell.cellNumX;
				int j      = cell.cellNumY;
				int center = i * n + j;
				int left   = (i - 1) * n + j;
				int right  = (i + 1) * n + j;
				int bottom = i * n + j - 1;
				int top    = i * n + j + 1;

				auto sx0 = gridCells[left].s;
				auto sx1 = gridCells[right].s;
				auto sy0 = gridCells[bottom].s;
				auto sy1 = gridCells[top].s;
				auto s   = sx0 + sx1 + sy0 + sy1;
				if (!isVeryCloseToZero(s)) {
					auto div = gridCells[right].u - gridCells[center].u + gridCells[top].v - gridCells[center].v;

					if (particleRestDensity > 0.0 && compensateDrift) {
						auto k           = 1.0;
						auto compression = cell.particleDensity - particleRestDensity;
						if (compression > 0.0)
							div = div - k * compression;
					}

					auto pValue = -div / s;
					pValue *= overRelaxation;

					gridCells[center].u -= sx0 * pValue;
					gridCells[right].u += sx1 * pValue;
					gridCells[center].v -= sy0 * pValue;
					gridCells[top].v += sy1 * pValue;
				}
			}
		};
		for_each(std::execution::seq, gridCells.begin(), gridCells.end(), parseIncompressibility);
	}
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
    bool compensateDrift,
    bool separateParticles,
    double obstacleRadius,
    bool instrument)
{
	std::vector<ExecutionLog> log;



	QElapsedTimer timer;

	if (instrument)
		timer.start();

	integrateParticles();
	if (instrument) {
		log.push_back(ExecutionLog { "integrateParticles", timer.nsecsElapsed() });
		timer.restart();
	}
	if (separateParticles) {
		pushParticlesApart(numParticleIters);
		if (instrument) {
			log.push_back(ExecutionLog { "pushParticlesApart", timer.nsecsElapsed() });
			timer.restart();
		}
	}
	handleParticleCollisions(obstacleX, obstacleY, obstacleRadius);
	if (instrument) {
		log.push_back(ExecutionLog { "handleParticleCollisions", timer.nsecsElapsed() });
		timer.restart();
	}
	transferVelocitiesToGrid();
	if (instrument) {
		log.push_back(ExecutionLog { "transferVelocitiesTrue", timer.nsecsElapsed() });
		timer.restart();
	}
	updateParticleDensity();
	if (instrument) {
		log.push_back(ExecutionLog { "updateParticleDensity", timer.nsecsElapsed() });
		timer.restart();
	}
	solveIncompressibility(numPressureIters, constants::dt, overRelaxation, compensateDrift);
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