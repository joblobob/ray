module;

#include <algorithm>

export module CellCalculations;

import Constants;

export struct Point {
	double x, y;
};

export struct Border {
	int minX = 0, minY = 0, maxX = constants::maxwidth, maxY = constants::maxwidth;
};

export constexpr int cellNumber(const Point& pt, const Border& border, double spacing)
{
	int xi = std::clamp((int)floor(pt.x * spacing), border.minX, border.maxX - 1);
	int yi = std::clamp((int)floor(pt.y * spacing), border.minY, border.maxY - 1);
	return xi * border.maxY + yi;
}