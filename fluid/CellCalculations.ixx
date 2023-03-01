module;

#include <algorithm>

export module CellCalculations;

import BaseStructures;
import Constants;

export constexpr int cellNumber(const double x, const double y, const Border& border, const double spacing)
{
	const int xi = std::clamp(int_floor(x * spacing), border.minX, border.maxX - 1);
	const int yi = std::clamp(int_floor(y * spacing), border.minY, border.maxY - 1);
	return xi * border.maxY + yi;
}

constexpr double my_abs(const double x) noexcept
{
	return x < 0.0 ? -x : x;
}

export constexpr bool isVeryCloseToZero(const double val)
{
	constexpr double epsilon = std::numeric_limits<double>::epsilon();
	return my_abs(val) <= epsilon * my_abs(val);
	// see Knuth section 4.2.2 pages 217-218
}