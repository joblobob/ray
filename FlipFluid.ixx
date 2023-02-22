module;

#include <algorithm>
#include <string>
#include <vector>

#include <QString>

export module FlipFluid;

import BaseStructures;
import Constants;
import CellCalculations;

export struct ExecutionLog {
	QString message;
	long long nsElapsed;
};

export struct FlipFluid {
	double density, fInvSpacing, particleRestDensity, pInvSpacing, particleRadius, h, obstacleVelX, obstacleVelY, obstacleX, obstacleY;
	int fNumX, fNumY, fNumCells, maxParticles, pNumX, pNumY, pNumCells;

	std::vector<Particle> particleMap;
	std::vector<Cell> gridCells;

	Border pBorder;
	Border fBorder;

	FlipFluid() = default;
	FlipFluid(double density, double width, double height, double spacing, double particleRadius, int maxParticles);

	void setupObstacle(double x, double y, bool reset);
	void transferVelocitiesToParticles();

	std::vector<ExecutionLog>
	simulate(double dt, double flipRatio, int numPressureIters, int numParticleIters, double overRelaxation, double obstacleRadius, bool instrument);

private:
	void parseVelocitiesParticles(Particle& particle);
};
