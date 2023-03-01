module;

#include <algorithm>
#include <execution>

export module ParticleDensity;

import BaseStructures;


export void updateParticleDensity(std::vector<Particle>& particleMap, std::vector<Cell>& gridCells)
{
	// reset density
	for_each(std::execution::seq, gridCells.begin(), gridCells.end(), [](Cell& c) { c.particleDensity = 0.0; });

	auto solveParticleDensity = [&](const Particle& p) {
		const double x = std::clamp(p.posX, constants::cellHeight, (constants::fNumX - 1) * constants::cellHeight);
		const double y = std::clamp(p.posY, constants::cellHeight, (constants::fNumY - 1) * constants::cellHeight);

		const int x0    = int_floor((x - constants::halfCellHeight) * constants::fInvSpacing);
		const double tx = ((x - constants::halfCellHeight) - x0 * constants::cellHeight) * constants::fInvSpacing;
		const int x1    = std::min(x0 + 1, constants::fNumX - 2);

		const int y0    = int_floor((y - constants::halfCellHeight) * constants::fInvSpacing);
		const double ty = ((y - constants::halfCellHeight) - y0 * constants::cellHeight) * constants::fInvSpacing;
		const int y1    = std::min(y0 + 1, constants::fNumY - 2);

		const double sx = 1.0 - tx;
		const double sy = 1.0 - ty;

		if (x0 < constants::fNumX && y0 < constants::fNumY)
			gridCells[x0 * constants::fNumY + y0].particleDensity += sx * sy;
		if (x1 < constants::fNumX && y0 < constants::fNumY)
			gridCells[x1 * constants::fNumY + y0].particleDensity += tx * sy;
		if (x1 < constants::fNumX && y1 < constants::fNumY)
			gridCells[x1 * constants::fNumY + y1].particleDensity += tx * ty;
		if (x0 < constants::fNumX && y1 < constants::fNumY)
			gridCells[x0 * constants::fNumY + y1].particleDensity += sx * ty;
	};
	for_each(std::execution::seq, particleMap.begin(), particleMap.end(), solveParticleDensity);
}

export void calculateParticleRestDensity(std::vector<Cell>& gridCells, double& particleRestDensity)
{
	double sum        = 0.0;
	int numFluidCells = 0;

	auto calculateSum = [&](Cell& cell) {
		if (cell.cellType == constants::CellType::Fluid) {
			sum += cell.particleDensity;
			numFluidCells++;
		}
	};
	for_each(std::execution::seq, gridCells.begin(), gridCells.end(), calculateSum);

	if (numFluidCells > 0)
		particleRestDensity = sum / numFluidCells;
}