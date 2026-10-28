module;

#include <algorithm>
#include <array>
#include <vector>
#include "ParallelForEach.h"

export module ParticleIncompressibility;

import BaseStructures;
import CellCalculations;
import Constants;

// constants
constexpr float overRelaxation { 1.9f };
constexpr int numPressureIters { 50 };
constexpr float k { 1.0f };

struct PressureCell {
	int index;
	int right;
	int top;
	float leftWeight;
	float rightWeight;
	float bottomWeight;
	float topWeight;
	float weightSum;
};

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
	std::array<std::vector<PressureCell>, 2> checkerboardCells;
	checkerboardCells[0].reserve(gridCells.size() / 2);
	checkerboardCells[1].reserve(gridCells.size() / 2);
	for (int i = 0; i < static_cast<int>(gridCells.size()); ++i) {
		if (gridCells[i].cellType == constants::CellType::Fluid) {
			fluidCellIndices.push_back(i);
			const Cell& cell = gridCells[i];
			const int x = cell.cellNumX;
			const int y = cell.cellNumY;
			const int left = (x - 1) * constants::fNumY + y;
			const int right = (x + 1) * constants::fNumY + y;
			const int bottom = x * constants::fNumY + y - 1;
			const int top = x * constants::fNumY + y + 1;
			const Cell& cellLeft = gridCells[left];
			const Cell& cellRight = gridCells[right];
			const Cell& cellBottom = gridCells[bottom];
			const Cell& cellTop = gridCells[top];
			const float leftWeight = cellLeft.s * inverseFaceDensity(cellLeft, cell);
			const float rightWeight = cellRight.s * inverseFaceDensity(cell, cellRight);
			const float bottomWeight = cellBottom.s * inverseFaceDensity(cellBottom, cell);
			const float topWeight = cellTop.s * inverseFaceDensity(cell, cellTop);
			checkerboardCells[(x + y) & 1].push_back({ i,
			    right,
			    top,
			    leftWeight,
			    rightWeight,
			    bottomWeight,
			    topWeight,
			    leftWeight + rightWeight + bottomWeight + topWeight });
		}
	}

	for (int iter = 0; iter < numPressureIters; iter++) {
		// Red/black Gauss-Seidel: cells of one color do not share face
		// velocities, so parallelize each color. The barrier between colors
		// preserves the pressure update dependency.
		for (const auto& colorCells : checkerboardCells) {
			fluidParallelForEach(colorCells.cbegin(), colorCells.cend(), [&](const PressureCell& pressureCell) {
				Cell& cell = gridCells[pressureCell.index];
				Cell& cellRight = gridCells[pressureCell.right];
				Cell& cellTop = gridCells[pressureCell.top];
				if (pressureCell.weightSum > 0.1f) {
					float div = cellRight.u - cell.u + cellTop.v - cell.v;

					const float compression = cell.particleDensity - particleRestDensity;
					if (compression > 0.0f)
						div = div - k * compression;

					const float pValue = (-div / pressureCell.weightSum) * overRelaxation;

					cell.u -= pressureCell.leftWeight * pValue;
					cellRight.u += pressureCell.rightWeight * pValue;
					cell.v -= pressureCell.bottomWeight * pValue;
					cellTop.v += pressureCell.topWeight * pValue;
				}
			});
		}
	}
}
