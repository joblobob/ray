module;

#include <algorithm>
#include <execution>

export module ParticleIncompressibility;

import BaseStructures;
import CellCalculations;


export void solveIncompressibility(std::vector<Cell>& gridCells,
    const int numIters,
    const int fNumY,
    const double density,
    const double h,
    const double dt,
    const double overRelaxation,
    const double particleRestDensity)
{
	const int n     = fNumY;
	const double cp = density * h / dt;

	for_each(std::execution::seq, gridCells.begin(), gridCells.end(), [&](Cell& cell) {
		cell.prevU = cell.u;
		cell.prevV = cell.v;
	});

	for (int iter = 0; iter < numIters; iter++) {
		auto parseIncompressibility = [&](Cell& cell) {
			if (cell.cellType == constants::CellType::Fluid) {
				const int i      = cell.cellNumX;
				const int j      = cell.cellNumY;
				const int center = i * n + j;
				const int left   = (i - 1) * n + j;
				const int right  = (i + 1) * n + j;
				const int bottom = i * n + j - 1;
				const int top    = i * n + j + 1;

				const double sx0 = gridCells[left].s;
				const double sx1 = gridCells[right].s;
				const double sy0 = gridCells[bottom].s;
				const double sy1 = gridCells[top].s;
				const double s   = sx0 + sx1 + sy0 + sy1;
				if (!isVeryCloseToZero(s)) {
					double div = gridCells[right].u - gridCells[center].u + gridCells[top].v - gridCells[center].v;

					const double k           = 1.0;
					const double compression = cell.particleDensity - particleRestDensity;
					if (compression > 0.0)
						div = div - k * compression;

					const double pValue = (-div / s) * overRelaxation;


					gridCells[center].u -= sx0 * pValue;
					gridCells[right].u += sx1 * pValue;
					gridCells[center].v -= sy0 * pValue;
					gridCells[top].v += sy1 * pValue;
				}
			}
		};
		for_each(std::execution::seq, gridCells.begin(), gridCells.end(), parseIncompressibility);
	}
}