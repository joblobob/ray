module;

#include <algorithm>
#include <execution>

export module CellColor;

import BaseStructures;
import CellCalculations;

void setSciColor(Cell& cell, double val, double minVal, double maxVal)
{
	val      = std::min(std::max(val, minVal), maxVal - 0.0001);
	double d = maxVal - minVal;
	val      = isVeryCloseToZero(d) ? 0.5 : (val - minVal) / d;
	double m = 0.25;
	int num  = floor(val / m);
	double s = (val - num * m) / m;
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
			auto d = cell.particleDensity;
			if (particleRestDensity > 0.0)
				d /= particleRestDensity;
			setSciColor(cell, d, 0.0, 2.0);
		}
	};
	for_each(std::execution::par_unseq, gridCells.begin(), gridCells.end(), calcCellColor);
}