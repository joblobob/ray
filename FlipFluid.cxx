module;

#include <QDebug>
#include <algorithm>
#include <execution>
#include <vector>

#include <QElapsedTimer.h>

module FlipFluid;

import CellColor;
import CellVelocities;
import ParticleIntegration;
import PushParticles;
import ParticleCollision;
import ParticleColor;
import ParticleDensity;
import ParticleIncompressibility;
import ParticleVelocities;

using namespace constants;

FlipFluid::FlipFluid(double density, double width, double height, double spacing, double particleRadius, int maxParticles) :


    // particles
    particleRestDensity(0.0), fBorder { .maxX = fNumX, .maxY = fNumY }, particleMap(maxParticles), gridCells(fNumCells), obstacle()
{
	h           = std::max(width / (double)fNumX, height / (double)fNumY);
	fInvSpacing = 1.0 / h;
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
	pushParticlesApart(particleMap, numParticleIters, pNumCells, pNumX, pNumY, maxParticles, pInvSpacing, r);
	if (instrument) {
		log.push_back(ExecutionLog { "pushParticlesApart", timer.nsecsElapsed() });
		timer.restart();
	}
	handleParticleCollisions(particleMap, fNumX, fNumY, fInvSpacing, obstacle, obstacleRadius, r);
	if (instrument) {
		log.push_back(ExecutionLog { "handleParticleCollisions", timer.nsecsElapsed() });
		timer.restart();
	}
	transferVelocitiesToGrid(particleMap, gridCells, h, fNumX, fNumY, fInvSpacing, fBorder);
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
	transferVelocitiesToParticles(particleMap, gridCells, h, fNumX, fNumY, fInvSpacing);
	if (instrument) {
		log.push_back(ExecutionLog { "transferVelocitiesFalse", timer.nsecsElapsed() });
		timer.restart();
	}


	updateParticleColors(particleMap, gridCells, fInvSpacing, fBorder, particleRestDensity);
	if (instrument) {
		log.push_back(ExecutionLog { "updateParticleColors", timer.nsecsElapsed() });
		timer.restart();
	}
	updateCellColors(gridCells, particleRestDensity);
	if (instrument) {
		log.push_back(ExecutionLog { "updateCellColors", timer.nsecsElapsed() });
		timer.restart();
	}
	return log;
}
