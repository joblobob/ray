module;

#include <QDebug>
#include <algorithm>
#include <execution>
#include <iostream>

export module ParticleDensity;

import BaseStructures;


export void updateParticleDensity(std::vector<Particle>& particleMap,
    std::vector<Cell>& gridCells,
    const int fNumX,
    const int fNumY,
    const double fInvSpacing,
    const double h)
{
	const int n     = fNumY;
	const double h1 = fInvSpacing;
	const double h2 = 0.5 * h;

	for_each(std::execution::seq, gridCells.begin(), gridCells.end(), [](Cell& c) { c.particleDensity = 0.0; });


	auto solveParticleDensity = [&](const Particle& p) {
		const double x = std::clamp(p.posX, h, (fNumX - 1) * h);
		const double y = std::clamp(p.posY, h, (fNumY - 1) * h);

		const int x0    = floor((x - h2) * h1);
		const double tx = ((x - h2) - x0 * h) * h1;
		const int x1    = std::min(x0 + 1, fNumX - 2);

		const int y0    = floor((y - h2) * h1);
		const double ty = ((y - h2) - y0 * h) * h1;
		const int y1    = std::min(y0 + 1, fNumY - 2);

		const double sx = 1.0 - tx;
		const double sy = 1.0 - ty;

		if (x0 < fNumX && y0 < fNumY)
			gridCells[x0 * n + y0].particleDensity += sx * sy;
		if (x1 < fNumX && y0 < fNumY)
			gridCells[x1 * n + y0].particleDensity += tx * sy;
		if (x1 < fNumX && y1 < fNumY)
			gridCells[x1 * n + y1].particleDensity += tx * ty;
		if (x0 < fNumX && y1 < fNumY)
			gridCells[x0 * n + y1].particleDensity += sx * ty;
	};
	for_each(std::execution::par_unseq, particleMap.begin(), particleMap.end(), solveParticleDensity);
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