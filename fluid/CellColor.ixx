module;

#include <algorithm>
#include <execution>

export module CellColor;

import BaseStructures;
import CellCalculations;

void setSciColor(Cell& cell, double val, const double minVal, const double maxVal)
{
	constexpr double epsilon = 0.0001;
	val                      = std::min(std::max(val, minVal), maxVal - epsilon);
	const double d           = maxVal - minVal;
	val                      = isVeryCloseToZero(d) ? 0.5 : (val - minVal) / d;
	constexpr double m       = 0.25;
	const unsigned int num { static_cast<unsigned int>(val / m) };
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

export void updateCellColors(std::vector<Cell>& gridCells, const double particleRestDensity)
{
	for_each(std::execution::seq, gridCells.begin(), gridCells.end(), [](Cell& c) {
		c.colorR = 0.0;
		c.colorG = 0.0;
		c.colorB = 0.0;
	});

	auto calcCellColor = [&](Cell& cell) {
		if (cell.cellType == constants::CellType::Solid) {
			cell.colorR = 0.5;
			cell.colorG = 0.5;
			cell.colorB = 0.5;
		} else if (cell.cellType == constants::CellType::Fluid) {
			double d = cell.particleDensity;
			if (particleRestDensity > 0.0)
				d /= particleRestDensity;
			setSciColor(cell, d, 0.0, 2.0);
		}
	};
	for_each(std::execution::par_unseq, gridCells.begin(), gridCells.end(), calcCellColor);
}