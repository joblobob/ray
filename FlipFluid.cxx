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


FlipFluid::FlipFluid(double width, double height, double spacing, double particleRadius, int maxParticles) :


    // particles
    particleRestDensity(0.0),
    fBorder { .maxX = constants::fNumX, .maxY = constants::fNumY },
    particleMap(maxParticles),
    gridCells(constants::fNumCells),
    obstacle()
{}


std::vector<ExecutionLog> FlipFluid::simulate(double dt, bool instrument)
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
	pushParticlesApart(particleMap);
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
	solveIncompressibility(gridCells, particleRestDensity);
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
