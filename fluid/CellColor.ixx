module;

#include <algorithm>
#include <execution>
#include <ranges>

export module CellColor;

import BaseStructures;
import CellCalculations;

void fluidColor(Cell& cell, float val, const float minVal, const float maxVal)
{
	constexpr float epsilon = 0.0001f;
	val                     = std::min(std::max(val, minVal), maxVal - epsilon);
	const float d           = maxVal - minVal;
	val                     = d < 0.1f ? 0.5f : (val - minVal) / d;
	constexpr float m       = 0.25f;
	const int num { static_cast<int>(val / m) };
	const float s = (val - num * m) / m;
	float r, g, b;

	switch (num) {
		case 0:
			r = 0.0f;
			g = s;
			b = 1.0f;
			break;
		case 1:
			r = 0.0f;
			g = 1.0f;
			b = 1.0f - s;
			break;
		case 2:
			r = s;
			g = 1.0f;
			b = 0.0f;
			break;
		case 3:
			r = 1.0f;
			g = 1.0f - s;
			b = 0.0f;
			break;
	}

	cell.colorR = r;
	cell.colorG = g;
	cell.colorB = b;
}

void solidCellColor(Cell& cell)
{
	cell.colorR = 0.5f;
	cell.colorG = 0.5f;
	cell.colorB = 0.5f;
}

void airCellColor(Cell& cell)
{
	cell.colorR = 0.0f;
	cell.colorG = 0.0f;
	cell.colorB = 0.0f;
}

export void updateCellColors(std::vector<Cell>& gridCells, const float particleRestDensity)
{
	// Solid
	//gridCells | std::views::filter(isSolid) | std::views::transform(solidCellColor) | std::ranges::to<std::vector>(); //slower, meh
	for (Cell& cell : gridCells | std::views::filter(isSolid)) {
		solidCellColor(cell);
	}
	// Air
	for (Cell& cell : gridCells | std::views::filter(isAir)) {
		airCellColor(cell);
	}

	// Fluid
	for (Cell& cell : gridCells | std::views::filter(isFluid)) {
		fluidColor(cell, cell.particleDensity / particleRestDensity, 0.0f, 2.0f);
	}
}