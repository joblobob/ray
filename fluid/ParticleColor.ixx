module;

#include <algorithm>
#include <execution>

export module ParticleColor;

import BaseStructures;
import CellCalculations;

export void updateParticleColors(std::vector<Particle>& particleMap,
    std::vector<Cell>& gridCells,
    const Border& fBorder,
    const float particleRestDensity)
{
	auto parseParticleColors = [&](Particle& particle) {
		const float s = 0.01f;

		particle.colorR = std::clamp(particle.colorR - s, 0.0f, 1.0f);
		particle.colorG = std::clamp(particle.colorG - s, 0.0f, 1.0f);
		particle.colorB = std::clamp(particle.colorB + s, 0.0f, 1.0f);

		const int cellNr = cellNumber(particle.posX, particle.posY, fBorder, constants::fInvSpacing);

		if (particleRestDensity > 0.0f) {
			const float relDensity = gridCells[cellNr].particleDensity / particleRestDensity;
			if (relDensity < 0.7f) {
				const float s   = 0.8f;
				particle.colorR = s;
				particle.colorG = s;
				particle.colorB = 1.0f;
			}
		}
	};
	for_each(std::execution::par_unseq, particleMap.begin(), particleMap.end(), parseParticleColors);
}