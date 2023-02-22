module;

#include <algorithm>
#include <execution>

export module ParticleColor;

import BaseStructures;
import CellCalculations;

export void updateParticleColors(std::vector<Particle>& particleMap,
    std::vector<Cell>& gridCells,
    const double fInvSpacing,
    const Border& fBorder,
    const double particleRestDensity)
{
	double h1 = fInvSpacing;

	auto parseParticleColors = [&](Particle& particle) {
		auto s = 0.01;

		particle.colorR = std::clamp(particle.colorR - s, 0.0, 1.0);
		particle.colorG = std::clamp(particle.colorG - s, 0.0, 1.0);
		particle.colorB = std::clamp(particle.colorB + s, 0.0, 1.0);

		int cellNr = cellNumber(particle.posX, particle.posY, fBorder, fInvSpacing);

		auto d0 = particleRestDensity;

		if (d0 > 0.0) {
			const auto relDensity = gridCells[cellNr].particleDensity / d0;
			if (relDensity < 0.7) {
				const auto s    = 0.8;
				particle.colorR = s;
				particle.colorG = s;
				particle.colorB = 1.0;
			}
		}
	};
	for_each(std::execution::par_unseq, particleMap.begin(), particleMap.end(), parseParticleColors);
}