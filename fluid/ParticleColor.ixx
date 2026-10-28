module;

#include <algorithm>
#include <vector>
#include "ParallelForEach.h"

export module ParticleColor;

import BaseStructures;
import CellCalculations;
import Constants;

export void updateParticleColors(ParticleStore& particleMap,
    std::vector<Cell>& gridCells,
    const Border& fBorder,
    const float particleRestDensity)
{
	auto parseParticleColors = [&](const std::size_t id) {
		const float s = 0.01f;

		particleMap.colorR[id] = std::clamp(particleMap.colorR[id] - s, 0.0f, 1.0f);
		particleMap.colorG[id] = std::clamp(particleMap.colorG[id] - s, 0.0f, 1.0f);
		particleMap.colorB[id] = std::clamp(particleMap.colorB[id] + s, 0.0f, 1.0f);

		const int cellNr = cellNumber(particleMap.posX[id], particleMap.posY[id], fBorder, constants::fInvSpacing);

		if (particleRestDensity > 0.0f) {
			const float relDensity = gridCells[cellNr].particleDensity / particleRestDensity;
			if (relDensity < 0.7f) {
				const float s   = 0.8f;
				particleMap.colorR[id] = s;
				particleMap.colorG[id] = s;
				particleMap.colorB[id] = 1.0f;
			}
		}
	};
	fluidParallelForEachIndex(particleMap.size(), parseParticleColors);
}
