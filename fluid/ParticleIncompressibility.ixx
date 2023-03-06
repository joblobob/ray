module;

#include <algorithm>
#include <execution>

export module ParticleIncompressibility;

import BaseStructures;
import CellCalculations;

// constants
constexpr float overRelaxation { 1.9f };
constexpr int numPressureIters { 50 };
constexpr float k { 1.0f };

// code
export void solveIncompressibility(std::vector<Cell>& gridCells, const float particleRestDensity)
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


				const float s = cellLeft.s + cellRight.s + cellBottom.s + cellTop.s;
				if (s > 0.1f) {
					float div = cellRight.u - cell.u + cellTop.v - cell.v;

					const float compression = cell.particleDensity - particleRestDensity;
					if (compression > 0.0f)
						div = div - k * compression;

					const float pValue = (-div / s) * overRelaxation;

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