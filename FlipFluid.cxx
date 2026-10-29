module;

#include <QDebug>
#include <algorithm>
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


FlipFluid::FlipFluid(float width, float height, float spacing, int maxParticles) :


    // particles
    particleRestDensity(0.0f),
    fBorder { .maxX = constants::fNumX, .maxY = constants::fNumY },
    particleMap(static_cast<std::size_t>(maxParticles)),
    gridCells(constants::fNumCells),
    obstacle()
{}


std::vector<ExecutionLog> FlipFluid::simulate(float dt, bool instrument, const int pressureIterations)
{
	std::vector<ExecutionLog> log;



	QElapsedTimer timer;

	if (instrument) {
		log.reserve(9);
		timer.start();
	}

	integrateParticles(particleMap, dt);
	if (instrument) {
		log.push_back(ExecutionLog { "integrateParticles", timer.nsecsElapsed() });
		timer.restart();
	}
	pushParticlesApart(particleMap, pushParticleBuffers);
	if (instrument) {
		log.push_back(ExecutionLog { "pushParticlesApart", timer.nsecsElapsed() });
		timer.restart();
	}
	handleParticleCollisions(particleMap, obstacle);
	if (instrument) {
		log.push_back(ExecutionLog { "handleParticleCollisions", timer.nsecsElapsed() });
		timer.restart();
	}
	transferVelocitiesToGrid(particleMap, gridCells, fBorder);
	if (instrument) {
		log.push_back(ExecutionLog { "transferVelocitiesTrue", timer.nsecsElapsed() });
		timer.restart();
	}
	updateParticleDensity(particleMap, gridCells);
	if (instrument) {
		log.push_back(ExecutionLog { "updateParticleDensity", timer.nsecsElapsed() });
		timer.restart();
	}
	if (particleRestDensity < 0.1) {
		calculateParticleRestDensity(gridCells, particleRestDensity);
	}
	solveIncompressibility(gridCells, particleRestDensity, fluidCellIndices, pressureIterations);
	if (instrument) {
		log.push_back(ExecutionLog { "solveIncompressibility", timer.nsecsElapsed() });
		timer.restart();
	}
	transferVelocitiesToParticles(particleMap, gridCells);
	if (instrument) {
		log.push_back(ExecutionLog { "transferVelocitiesFalse", timer.nsecsElapsed() });
		timer.restart();
	}


	updateParticleColors(particleMap, gridCells, fBorder, particleRestDensity);
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
