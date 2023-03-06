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
				const int i = cell.cellNumX;
				const int j = cell.cellNumY;

				Cell& cellLeft   = gridCells[(i - 1) * constants::fNumY + j]; //left
				Cell& cellRight  = gridCells[(i + 1) * constants::fNumY + j]; //right
				Cell& cellBottom = gridCells[i * constants::fNumY + j - 1];   //bottom
				Cell& cellTop    = gridCells[i * constants::fNumY + j + 1];   //top


				const double s = cellLeft.s + cellRight.s + cellBottom.s + cellTop.s;
				if (s > 0.1) {
					double div = cellRight.u - cell.u + cellTop.v - cell.v;

					const double compression = cell.particleDensity - particleRestDensity;
					if (compression > 0.0)
						div = div - k * compression;

					const double pValue = (-div / s) * overRelaxation;

					cell.u -= cellLeft.s * pValue;
					cellRight.u += cellRight.s * pValue;
					cell.v -= cellBottom.s * pValue;
					cellTop.v += cellTop.s * pValue;
				}
			}
		};
		for_each(std::execution::seq, gridCells.begin(), gridCells.end(), parseIncompressibility);
	}
}