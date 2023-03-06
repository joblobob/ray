module;

#include <algorithm>

export module CellCalculations;

import BaseStructures;
import Constants;

export constexpr int cellNumber(const double x, const double y, const Border& border, const double spacing)
{
	const int xi = std::clamp(static_cast<int>(x * spacing), border.minX, border.maxX - 1);
	const int yi = std::clamp(static_cast<int>(y * spacing), border.minY, border.maxY - 1);
	return xi * border.maxY + yi;
}

constexpr double my_abs(const double x) noexcept
{
	return x < 0.0 ? -x : x;
}
