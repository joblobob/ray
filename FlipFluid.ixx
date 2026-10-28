module;

#include <algorithm>
#include <string>
#include <vector>

export module FlipFluid;

import BaseStructures;
import Constants;
import CellCalculations;
import Obstacle;
import PushParticles;

export struct ExecutionLog {
	std::string message;
	long long nsElapsed;
};

export struct FlipFluid {
	float particleRestDensity;

	Obstacle obstacle;
	ParticleStore particleMap;
	std::vector<Cell> gridCells;
	std::vector<int> fluidCellIndices;
	PushParticlesBuffers pushParticleBuffers;

	Border fBorder;

	FlipFluid() = default;
	FlipFluid(float width, float height, float spacing, int maxParticles);

	std::vector<ExecutionLog> simulate(float dt, bool instrument);
};
