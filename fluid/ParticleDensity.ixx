module;

#include <algorithm>
#include <vector>

export module ParticleDensity;

import BaseStructures;
import Constants;


export void updateParticleDensity(const ParticleStore& particleMap, std::vector<Cell>& gridCells)
{
	// reset density
	std::for_each(gridCells.begin(), gridCells.end(), [](Cell& c) {
		c.particleDensity = 0.0f;
		c.particleMassDensity = 0.0f;
	});

	auto solveParticleDensity = [&](const int id) {
		const float mass = constants::particleDensityRatioForId(id);
		const float x = std::clamp(particleMap.posX[id], constants::cellHeight, (constants::fNumX - 1) * constants::cellHeight);
		const float y = std::clamp(particleMap.posY[id], constants::cellHeight, (constants::fNumY - 1) * constants::cellHeight);

		const int x0   = static_cast<int>((x - constants::halfCellHeight) * constants::fInvSpacing);
		const float tx = ((x - constants::halfCellHeight) - x0 * constants::cellHeight) * constants::fInvSpacing;
		const int x1   = std::min(x0 + 1, constants::fNumX - 2);

		const int y0   = static_cast<int>((y - constants::halfCellHeight) * constants::fInvSpacing);
		const float ty = ((y - constants::halfCellHeight) - y0 * constants::cellHeight) * constants::fInvSpacing;
		const int y1   = std::min(y0 + 1, constants::fNumY - 2);

		const float sx = 1.0 - tx;
		const float sy = 1.0 - ty;

		if (x0 < constants::fNumX && y0 < constants::fNumY)
		{
			Cell& cell = gridCells[x0 * constants::fNumY + y0];
			const float weight = sx * sy;
			cell.particleDensity += weight;
			cell.particleMassDensity += weight * mass;
		}
		if (x1 < constants::fNumX && y0 < constants::fNumY)
		{
			Cell& cell = gridCells[x1 * constants::fNumY + y0];
			const float weight = tx * sy;
			cell.particleDensity += weight;
			cell.particleMassDensity += weight * mass;
		}
		if (x1 < constants::fNumX && y1 < constants::fNumY)
		{
			Cell& cell = gridCells[x1 * constants::fNumY + y1];
			const float weight = tx * ty;
			cell.particleDensity += weight;
			cell.particleMassDensity += weight * mass;
		}
		if (x0 < constants::fNumX && y1 < constants::fNumY)
		{
			Cell& cell = gridCells[x0 * constants::fNumY + y1];
			const float weight = sx * ty;
			cell.particleDensity += weight;
			cell.particleMassDensity += weight * mass;
		}
	};
	for (int id = 0; id < static_cast<int>(particleMap.size()); ++id)
		solveParticleDensity(id);
}

export void calculateParticleRestDensity(std::vector<Cell>& gridCells, float& particleRestDensity)
{
	float sum         = 0.0f;
	int numFluidCells = 0;

	auto calculateSum = [&](Cell& cell) {
		if (cell.cellType == constants::CellType::Fluid) {
			sum += cell.particleDensity;
			numFluidCells++;
		}
	};
	std::for_each(gridCells.begin(), gridCells.end(), calculateSum);

	if (numFluidCells > 0)
		particleRestDensity = sum / numFluidCells;
}
