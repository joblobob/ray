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
				const int left   = (i - 1) * constants::fNumY + j;
				const int right  = (i + 1) * constants::fNumY + j;
				const int bottom = i * constants::fNumY + j - 1;
				const int top    = i * constants::fNumY + j + 1;

				auto cellLeft   = gridCells[left];
				auto& cellRight = gridCells[right];
				auto cellBottom = gridCells[bottom];
				auto& cellTop   = gridCells[top];


				const double sx0 = cellLeft.s;
				const double sx1 = cellRight.s;
				const double sy0 = cellBottom.s;
				const double sy1 = cellTop.s;
				const double s   = sx0 + sx1 + sy0 + sy1;
				if (!isVeryCloseToZero(s)) {
					double div = cellRight.u - cell.u + cellTop.v - cell.v;

					const double compression = cell.particleDensity - particleRestDensity;
					if (compression > 0.0)
						div = div - k * compression;

					const double pValue = (-div / s) * overRelaxation;

					cell.u -= sx0 * pValue;
					cellRight.u += sx1 * pValue;
					cell.v -= sy0 * pValue;
					cellTop.v += sy1 * pValue;
				}
			}
		};
		for_each(std::execution::seq, gridCells.begin(), gridCells.end(), parseIncompressibility);
	}
}