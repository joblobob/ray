module;

#include <algorithm>
#include <string>
#include <vector>

#include <QString>

export module FlipFluid;

import BaseStructures;
import Constants;
import CellCalculations;
import Obstacle;

export struct ExecutionLog {
	QString message;
	long long nsElapsed;
};

export struct FlipFluid {
	double fInvSpacing, particleRestDensity, h;

	Obstacle obstacle;
	std::vector<Particle> particleMap;
	std::vector<Cell> gridCells;

	Border fBorder;

	FlipFluid() = default;
	FlipFluid(double density, double width, double height, double spacing, double particleRadius, int maxParticles);

	std::vector<ExecutionLog>
	simulate(double dt, double flipRatio, int numPressureIters, int numParticleIters, double overRelaxation, double obstacleRadius, bool instrument);
};
