module;

#include <algorithm>
#include <execution>

export module ParticleIncompressibility;

import BaseStructures;
import CellCalculations;

// constants
constexpr double overRelaxation { 1.9 };
constexpr double numPressureIters { 50 };
constexpr double k { 1.0 };

// code
export void solveIncompressibility(std::vector<Cell>& gridCells, const double particleRestDensity)
{
	for_each(std::execution::seq, gridCells.begin(), gridCells.end(), [&](Cell& cell) {
		cell.prevU = cell.u;
		cell.prevV = cell.v;
	});

	for (int iter = 0; iter < numPressureIters; iter++) {
		auto parseIncompressibility = [&](Cell& cell) {
			if (cell.cellType == constants::CellType::Fluid) {
				const int i      = cell.cellNumX;
				const int j      = cell.cellNumY;
				const int center = i * constants::fNumY + j;
				const int left   = (i - 1) * constants::fNumY + j;
				const int right  = (i + 1) * constants::fNumY + j;
				const int bottom = i * constants::fNumY + j - 1;
				const int top    = i * constants::fNumY + j + 1;

				const double sx0 = gridCells[left].s;
				const double sx1 = gridCells[right].s;
				const double sy0 = gridCells[bottom].s;
				const double sy1 = gridCells[top].s;
				const double s   = sx0 + sx1 + sy0 + sy1;
				if (!isVeryCloseToZero(s)) {
					double div = gridCells[right].u - gridCells[center].u + gridCells[top].v - gridCells[center].v;

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