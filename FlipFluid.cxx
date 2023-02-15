module;

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
    u(fNumCells),
    // fluid
    v(fNumCells),
    du(fNumCells),
    dv(fNumCells),
    prevU(fNumCells),
    prevV(fNumCells),
    s(fNumCells),
    cellType(fNumCells),
    cellColor(3 * fNumCells),
    particleRadius(particleRadius),
    // particles
    maxParticles(maxParticles),
    particleDensity(fNumCells),
    particleRestDensity(0.0),
    pInvSpacing(1.0 / (2.2 * particleRadius)),
    pNumX(floor(width * pInvSpacing) + 1),
    pNumY(floor(height * pInvSpacing) + 1),
    pNumCells(pNumX * pNumY),
    numCellParticles(pNumCells),
    firstCellParticle(pNumCells + 1),
    cellParticleIds(maxParticles),
    pBorder { .maxX = pNumX, .maxY = pNumY },
    fBorder { .maxX = fNumX, .maxY = fNumY }
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
			s[i * n + j] = 1.0;

			auto dx = (i + 0.5) * h - x;
			auto dy = (j + 0.5) * h - y;

			if (dx * dx + dy * dy < r * r) {
				s[i * n + j]       = 0.0;
				u[i * n + j]       = vx;
				u[(i + 1) * n + j] = vx;
				v[i * n + j]       = vy;
				v[i * n + j + 1]   = vy;
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
		particle.particleVelY += constants::integ; //dt * gravity
		particle.particlePosX += particle.particleVelX * constants::dt;
		particle.particlePosY += particle.particleVelY * constants::dt;
	};
	std::for_each(std::execution::seq, particleMap.begin(), particleMap.end(), integration);
}

void FlipFluid::pushParticlesApart(int numIters)
{
	std::fill(numCellParticles.begin(), numCellParticles.end(), 0);

	// count particles per cell
	auto countParticlesInCell = [&](const Particle& particle) {
		numCellParticles[cellNumber(particle.particlePosX, particle.particlePosY, pBorder, pInvSpacing)]++;
	};
	for_each(std::execution::par_unseq, particleMap.begin(), particleMap.end(), countParticlesInCell);

	// partial sums
	auto first = 0;

	for (int i = 0; i < pNumCells; i++) {
		first += numCellParticles[i];
		firstCellParticle[i] = first;
	}
	firstCellParticle[pNumCells] = first; // guard

	// fill particles into cells
	auto fillParticleIntoCell = [&](const Particle& p) {
		const int cellNr = cellNumber(p.particlePosX, p.particlePosY, pBorder, pInvSpacing);
		firstCellParticle[cellNr]--;
		cellParticleIds[firstCellParticle[cellNr]] = p.particleId;
	};
	for_each(std::execution::par_unseq, particleMap.begin(), particleMap.end(), fillParticleIntoCell);


	// push particles apart

	const auto minDist  = 2.0 * particleRadius;
	const auto minDist2 = minDist * minDist;

	auto pushParticles = [&](auto& particle) {
		const auto px = particle.particlePosX;
		const auto py = particle.particlePosY;

		const int pxi = floor(px * pInvSpacing);
		const int pyi = floor(py * pInvSpacing);
		const int x0  = std::max(pxi - 1, 0);
		const int y0  = std::max(pyi - 1, 0);
		const int x1  = std::min(pxi + 1, pNumX - 1);
		const int y1  = std::min(pyi + 1, pNumY - 1);

		for (int xi = x0; xi <= x1; xi++) {
			for (int yi = y0; yi <= y1; yi++) {
				const int cellNr = xi * pNumY + yi;
				const int first  = firstCellParticle[cellNr];
				const int last   = firstCellParticle[cellNr + 1];
				for (int j = first; j < last; j++) {
					int id = cellParticleIds[j];
					if (id == particle.particleId)
						continue;
					auto& particleAtId { particleMap[id] };
					const auto qx = particleAtId.particlePosX;
					const auto qy = particleAtId.particlePosY;

					auto dx       = qx - px;
					auto dy       = qy - py;
					const auto d2 = dx * dx + dy * dy;
					if (d2 > minDist2 || isVeryCloseToZero(d2))
						continue;
					const auto d = sqrt(d2);
					const auto s = 0.5 * (minDist - d) / d;
					dx *= s;
					dy *= s;
					particle.particlePosX -= dx;
					particle.particlePosY -= dy;
					particleAtId.particlePosX += dx;
					particleAtId.particlePosY += dy;

					// diffuse colors

					calcColor(particle.particleColorR, particleAtId.particleColorR);
					calcColor(particle.particleColorG, particleAtId.particleColorG);
					calcColor(particle.particleColorB, particleAtId.particleColorB);
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
		auto x = particle.particlePosX;
		auto y = particle.particlePosY;

		const auto dx = x - obstacleX;
		const auto dy = y - obstacleY;
		const auto d2 = dx * dx + dy * dy;

		// obstacle collision

		if (d2 < minDist2) {
			particle.particleVelX = obstacleVelX;
			particle.particleVelY = obstacleVelY;
		}
		// wall collisions

		if (x < minX) {
			x                     = minX;
			particle.particleVelX = 0.0;
		} else if (x > maxX) {
			x                     = maxX;
			particle.particleVelX = 0.0;
		}
		if (y < minY) {
			y                     = minY;
			particle.particleVelY = 0.0;
		} else if (y > maxY) {
			y                     = maxY;
			particle.particleVelY = 0.0;
		}
		particle.particlePosX = x;
		particle.particlePosY = y;
	};
	for_each(std::execution::seq, particleMap.begin(), particleMap.end(), particleCollision);
}

void FlipFluid::updateParticleDensity()
{
	const int n     = fNumY;
	const auto h1   = fInvSpacing;
	const double h2 = 0.5 * h;

	std::fill(particleDensity.begin(), particleDensity.end(), 0.0);

	auto solveParticleDensity = [&](const Particle& p) {
		const auto x = std::clamp(p.particlePosX, h, (fNumX - 1) * h);
		const auto y = std::clamp(p.particlePosY, h, (fNumY - 1) * h);

		const int x0  = floor((x - h2) * h1);
		const auto tx = ((x - h2) - x0 * h) * h1;
		const int x1  = std::min(x0 + 1, fNumX - 2);

		const int y0  = floor((y - h2) * h1);
		const auto ty = ((y - h2) - y0 * h) * h1;
		const int y1  = std::min(y0 + 1, fNumY - 2);

		const auto sx = 1.0 - tx;
		const auto sy = 1.0 - ty;

		if (x0 < fNumX && y0 < fNumY)
			particleDensity[x0 * n + y0] += sx * sy;
		if (x1 < fNumX && y0 < fNumY)
			particleDensity[x1 * n + y0] += tx * sy;
		if (x1 < fNumX && y1 < fNumY)
			particleDensity[x1 * n + y1] += tx * ty;
		if (x0 < fNumX && y1 < fNumY)
			particleDensity[x0 * n + y1] += sx * ty;
	};
	for_each(std::execution::par_unseq, particleMap.begin(), particleMap.end(), solveParticleDensity);


	if (isVeryCloseToZero(particleRestDensity)) {
		auto sum           = 0.0;
		auto numFluidCells = 0;

		for (auto i = 0; i < fNumCells; i++) {
			if (cellType[i] == constants::CellType::Fluid) {
				sum += particleDensity[i];
				numFluidCells++;
			}
		}

		if (numFluidCells > 0)
			particleRestDensity = sum / numFluidCells;
	}
}

void FlipFluid::transferVelocities(bool toGrid, double flipRatio)
{
	auto n    = fNumY;
	auto h1   = fInvSpacing;
	double h2 = 0.5 * h;

	if (toGrid) {
		prevU = u;
		prevV = v;

		std::fill(du.begin(), du.end(), 0.0);
		std::fill(dv.begin(), dv.end(), 0.0);
		std::fill(u.begin(), u.end(), 0.0);
		std::fill(v.begin(), v.end(), 0.0);

		for (int i = 0; i < fNumCells; i++)
			cellType[i] = isVeryCloseToZero(s[i]) ? constants::CellType::Solid : constants::CellType::Air;

		auto setCellType = [&](const Particle& p) {
			int cellNr = cellNumber({ p.particlePosX, p.particlePosY }, fBorder, fInvSpacing);
			if (cellType[cellNr] == constants::CellType::Air)
				cellType[cellNr] = constants::CellType::Fluid;
		};
		for_each(std::execution::par_unseq, particleMap.begin(), particleMap.end(), setCellType);
	}


	for (int component = 0; component < 2; component++) {
		auto dx = component == 0 ? 0.0 : h2;
		auto dy = component == 0 ? h2 : 0.0;

		std::vector<double>& f     = component == 0 ? u : v;
		std::vector<double>& prevF = component == 0 ? prevU : prevV;
		std::vector<double>& d     = component == 0 ? du : dv;

		auto parseVelocities = [&](Particle& particle) {
			const auto x = std::clamp(particle.particlePosX, h, (fNumX - 1) * h);
			const auto y = std::clamp(particle.particlePosY, h, (fNumY - 1) * h);

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

			if (toGrid) {
				auto pv = component == 0 ? particle.particleVelX : particle.particleVelY;
				f[nr0] += pv * d0;
				d[nr0] += d0;
				f[nr1] += pv * d1;
				d[nr1] += d1;
				f[nr2] += pv * d2;
				d[nr2] += d2;
				f[nr3] += pv * d3;
				d[nr3] += d3;
			} else {
				int offset  = component == 0 ? n : 1;
				auto valid0 = cellType[nr0] != constants::CellType::Air || cellType[nr0 - offset] != constants::CellType::Air ? 1.0 : 0.0;
				auto valid1 = cellType[nr1] != constants::CellType::Air || cellType[nr1 - offset] != constants::CellType::Air ? 1.0 : 0.0;
				auto valid2 = cellType[nr2] != constants::CellType::Air || cellType[nr2 - offset] != constants::CellType::Air ? 1.0 : 0.0;
				auto valid3 = cellType[nr3] != constants::CellType::Air || cellType[nr3 - offset] != constants::CellType::Air ? 1.0 : 0.0;

				auto v = component == 0 ? particle.particleVelX : particle.particleVelY;
				auto d = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

				if (d > 0.0) {
					auto picV = (valid0 * d0 * f[nr0] + valid1 * d1 * f[nr1] + valid2 * d2 * f[nr2] + valid3 * d3 * f[nr3]) / d;
					auto corr = (valid0 * d0 * (f[nr0] - prevF[nr0]) + valid1 * d1 * (f[nr1] - prevF[nr1]) + valid2 * d2 * (f[nr2] - prevF[nr2]) +
					                valid3 * d3 * (f[nr3] - prevF[nr3])) /
					            d;
					auto flipV = v + corr;

					auto newVelValue = (1.0 - flipRatio) * picV + flipRatio * flipV;
					if (component == 0)
						particle.particleVelX = newVelValue;
					else
						particle.particleVelY = newVelValue;
				}
			}
		};

		for_each(std::execution::par_unseq, particleMap.begin(), particleMap.end(), parseVelocities);


		if (toGrid) {
			for (int i = 0; i < f.size(); i++) {
				if (d[i] > 0.0)
					f[i] /= d[i];
			}

			// restore solid cells

			for (int i = 0; i < fNumX; i++) {
				for (int j = 0; j < fNumY; j++) {
					auto solid = cellType[i * n + j] == constants::CellType::Solid;
					if (solid || (i > 0 && cellType[(i - 1) * n + j] == constants::CellType::Solid))
						u[i * n + j] = prevU[i * n + j];
					if (solid || (j > 0 && cellType[i * n + j - 1] == constants::CellType::Solid))
						v[i * n + j] = prevV[i * n + j];
				}
			}
		}
	}
}

void FlipFluid::solveIncompressibility(int numIters, double dt, double overRelaxation, bool compensateDrift)
{
	prevU = u;
	prevV = v;

	auto n  = fNumY;
	auto cp = density * h / dt;


	for (int iter = 0; iter < numIters; iter++) {
		for (int i = 1; i < fNumX - 1; i++) {
			for (int j = 1; j < fNumY - 1; j++) {
				if (cellType[i * n + j] != constants::CellType::Fluid)
					continue;

				int center = i * n + j;
				int left   = (i - 1) * n + j;
				int right  = (i + 1) * n + j;
				int bottom = i * n + j - 1;
				int top    = i * n + j + 1;

				auto sx0 = s[left];
				auto sx1 = s[right];
				auto sy0 = s[bottom];
				auto sy1 = s[top];
				auto s   = sx0 + sx1 + sy0 + sy1;
				if (isVeryCloseToZero(s))
					continue;

				auto div = u[right] - u[center] + v[top] - v[center];

				if (particleRestDensity > 0.0 && compensateDrift) {
					auto k           = 1.0;
					auto compression = particleDensity[i * n + j] - particleRestDensity;
					if (compression > 0.0)
						div = div - k * compression;
				}

				auto pValue = -div / s;
				pValue *= overRelaxation;

				u[center] -= sx0 * pValue;
				u[right] += sx1 * pValue;
				v[center] -= sy0 * pValue;
				v[top] += sy1 * pValue;
			}
		}
	}
}

void FlipFluid::updateParticleColors()
{
	auto h1 = fInvSpacing;

	auto parseParticleColors = [&](Particle& particle) {
		auto s = 0.01;

		particle.particleColorR = std::clamp(particle.particleColorR - s, 0.0, 1.0);
		particle.particleColorG = std::clamp(particle.particleColorG - s, 0.0, 1.0);
		particle.particleColorB = std::clamp(particle.particleColorB + s, 0.0, 1.0);

		int cellNr = cellNumber({ particle.particlePosX, particle.particlePosY }, fBorder, fInvSpacing);

		auto d0 = particleRestDensity;

		if (d0 > 0.0) {
			const auto relDensity = particleDensity[cellNr] / d0;
			if (relDensity < 0.7) {
				const auto s            = 0.8;
				particle.particleColorR = s;
				particle.particleColorG = s;
				particle.particleColorB = 1.0;
			}
		}
	};
	for_each(std::execution::par_unseq, particleMap.begin(), particleMap.end(), parseParticleColors);
}

void FlipFluid::setSciColor(int cellNr, double val, double minVal, double maxVal)
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

	cellColor[3 * cellNr]     = r;
	cellColor[3 * cellNr + 1] = g;
	cellColor[3 * cellNr + 2] = b;
}

void FlipFluid::updateCellColors()
{
	std::fill(cellColor.begin(), cellColor.end(), 0.0);

	for (int i = 0; i < fNumCells; i++) {
		if (cellType[i] == constants::CellType::Solid) {
			cellColor[3 * i]     = 0.5;
			cellColor[3 * i + 1] = 0.5;
			cellColor[3 * i + 2] = 0.5;
		} else if (cellType[i] == constants::CellType::Fluid) {
			auto d = particleDensity[i];
			if (particleRestDensity > 0.0)
				d /= particleRestDensity;
			setSciColor(i, d, 0.0, 2.0);
		}
	}
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
	transferVelocities(true, flipRatio);
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
	transferVelocities(false, flipRatio);
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
