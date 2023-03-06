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

constexpr float my_abs(const float x) noexcept
{
	return x < 0.0f ? -x : x;
}
