module;

#include <algorithm>

export module CellCalculations;

import BaseStructures;
import Constants;


export constexpr int cellNumber(const double& x, const double& y, const Border& border, double spacing)
{
	int xi = std::clamp((int)floor(x * spacing), border.minX, border.maxX - 1);
	int yi = std::clamp((int)floor(y * spacing), border.minY, border.maxY - 1);
	return xi * border.maxY + yi;
}

export inline bool isVeryCloseToZero(double val)
{
	constexpr double epsilon = std::numeric_limits<double>::epsilon();
	return std::abs(val) <= epsilon * std::abs(val);
	// see Knuth section 4.2.2 pages 217-218
}