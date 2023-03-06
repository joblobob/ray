module;

#include <algorithm>
#include <execution>
#include <ranges>

export module CellColor;

import BaseStructures;
import CellCalculations;

void fluidColor(Cell& cell, double val, const double minVal, const double maxVal)
{
	constexpr double epsilon = 0.0001;
	val                      = std::min(std::max(val, minVal), maxVal - epsilon);
	const double d           = maxVal - minVal;
	val                      = d < 0.1 ? 0.5 : (val - minVal) / d;
	constexpr double m       = 0.25;
	const int num { static_cast<int>(val / m) };
	const double s = (val - num * m) / m;
	double r, g, b;

	switch (num) {
		case 0:
			r = 0.0;
			g = s;
			b = 1.0;
			break;
		case 1:
			r = 0.0;
			g = 1.0;
			b = 1.0 - s;
			break;
		case 2:
			r = s;
			g = 1.0;
			b = 0.0;
			break;
		case 3:
			r = 1.0;
			g = 1.0 - s;
			b = 0.0;
			break;
	}

	cell.colorR = r;
	cell.colorG = g;
	cell.colorB = b;
}

bool isSolid(const Cell& cell)
{
	return cell.cellType == constants::CellType::Solid;
}
bool isFluid(const Cell& cell)
{
	return cell.cellType == constants::CellType::Fluid;
}
bool isAir(const Cell& cell)
{
	return cell.cellType == constants::CellType::Air;
}

void solidCellColor(Cell& cell)
{
	cell.colorR = 0.5;
	cell.colorG = 0.5;
	cell.colorB = 0.5;
}

void airCellColor(Cell& cell)
{
	cell.colorR = 0.0;
	cell.colorG = 0.0;
	cell.colorB = 0.0;
}

export void updateCellColors(std::vector<Cell>& gridCells, const double particleRestDensity)
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
		fluidColor(cell, cell.particleDensity / particleRestDensity, 0.0, 2.0);
	}
}