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
	float particleRestDensity;

	Obstacle obstacle;
	std::vector<Particle> particleMap;
	std::vector<Cell> gridCells;

	Border fBorder;

	FlipFluid() = default;
	FlipFluid(float width, float height, float spacing, float particleRadius, int maxParticles);

	std::vector<ExecutionLog> simulate(float dt, bool instrument);
};
