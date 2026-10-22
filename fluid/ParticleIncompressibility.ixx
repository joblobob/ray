module;

#include <algorithm>
#include <vector>

export module ParticleIncompressibility;

import BaseStructures;
import CellCalculations;
import Constants;

// constants
constexpr float overRelaxation { 1.9f };
constexpr int numPressureIters { 50 };
constexpr float k { 1.0f };

constexpr float cellMaterialDensity(const Cell& cell)
{
	return cell.particleDensity > 1.0e-6f ? cell.particleMassDensity / cell.particleDensity : 1.0f;
}

constexpr float inverseFaceDensity(const Cell& first, const Cell& second)
{
	const bool firstHasParticles = first.particleDensity > 1.0e-6f;
	const bool secondHasParticles = second.particleDensity > 1.0e-6f;
	if (firstHasParticles && secondHasParticles)
		return 2.0f / (cellMaterialDensity(first) + cellMaterialDensity(second));
	if (firstHasParticles)
		return 1.0f / cellMaterialDensity(first);
	if (secondHasParticles)
		return 1.0f / cellMaterialDensity(second);
	return 1.0f;
}

// code
export void solveIncompressibility(std::vector<Cell>& gridCells, const float particleRestDensity, std::vector<int>& fluidCellIndices)
{
	std::for_each(gridCells.begin(), gridCells.end(), [&](Cell& cell) {
		cell.prevU = cell.u;
		cell.prevV = cell.v;
	});
	fluidCellIndices.clear();
	for (int i = 0; i < static_cast<int>(gridCells.size()); ++i) {
		if (gridCells[i].cellType == constants::CellType::Fluid)
			fluidCellIndices.push_back(i);
	}

	for (int iter = 0; iter < numPressureIters; iter++) {
		for (const int cellIndex : fluidCellIndices) {
			Cell& cell = gridCells[cellIndex];
			const int i = cell.cellNumX;
			const int j = cell.cellNumY;

			Cell& cellLeft   = gridCells[(i - 1) * constants::fNumY + j]; //left
			Cell& cellRight  = gridCells[(i + 1) * constants::fNumY + j]; //right
			Cell& cellBottom = gridCells[i * constants::fNumY + j - 1];   //bottom
			Cell& cellTop    = gridCells[i * constants::fNumY + j + 1];   //top


			const float leftWeight = cellLeft.s * inverseFaceDensity(cellLeft, cell);
			const float rightWeight = cellRight.s * inverseFaceDensity(cell, cellRight);
			const float bottomWeight = cellBottom.s * inverseFaceDensity(cellBottom, cell);
			const float topWeight = cellTop.s * inverseFaceDensity(cell, cellTop);
			const float weightSum = leftWeight + rightWeight + bottomWeight + topWeight;
			if (weightSum > 0.1f) {
				float div = cellRight.u - cell.u + cellTop.v - cell.v;

				const float compression = cell.particleDensity - particleRestDensity;
				if (compression > 0.0f)
					div = div - k * compression;

				const float pValue = (-div / weightSum) * overRelaxation;

				cell.u -= leftWeight * pValue;
				cellRight.u += rightWeight * pValue;
				cell.v -= bottomWeight * pValue;
				cellTop.v += topWeight * pValue;
			}
		}
	}
}
