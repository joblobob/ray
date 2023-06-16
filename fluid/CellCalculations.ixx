module;

#include <algorithm>

export module CellCalculations;

import BaseStructures;
import Constants;

export constexpr int cellNumber(const float x, const float y, const Border& border, const float spacing)
{
	const int xi = std::clamp(static_cast<int>(x * spacing), border.minX, border.maxX - 1);
	const int yi = std::clamp(static_cast<int>(y * spacing), border.minY, border.maxY - 1);
	return xi * border.maxY + yi;
}


export bool isSolid(const Cell& cell)
{
	return cell.cellType == constants::CellType::Solid;
}
export bool isFluid(const Cell& cell)
{
	return cell.cellType == constants::CellType::Fluid;
}
export bool isAir(const Cell& cell)
{
	return cell.cellType == constants::CellType::Air;
}
