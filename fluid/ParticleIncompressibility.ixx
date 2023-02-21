module;

#include <QDebug>
#include <algorithm>
#include <execution>
#include <iostream>

export module ParticleIncompressibility;

import BaseStructures;
import CellCalculations;


export void solveIncompressibility(std::vector<Cell>& gridCells,
    int numIters,
    int fNumY,
    double density,
    double h,
    double dt,
    double overRelaxation,
    double particleRestDensity)
{
	auto n  = fNumY;
	auto cp = density * h / dt;

	for_each(std::execution::seq, gridCells.begin(), gridCells.end(), [&](Cell& cell) {
		cell.prevU = cell.u;
		cell.prevV = cell.v;
	});

	for (int iter = 0; iter < numIters; iter++) {
		auto parseIncompressibility = [&](Cell& cell) {
			if (cell.cellType == constants::CellType::Fluid) {
				int i      = cell.cellNumX;
				int j      = cell.cellNumY;
				int center = i * n + j;
				int left   = (i - 1) * n + j;
				int right  = (i + 1) * n + j;
				int bottom = i * n + j - 1;
				int top    = i * n + j + 1;

				auto sx0 = gridCells[left].s;
				auto sx1 = gridCells[right].s;
				auto sy0 = gridCells[bottom].s;
				auto sy1 = gridCells[top].s;
				auto s   = sx0 + sx1 + sy0 + sy1;
				if (!isVeryCloseToZero(s)) {
					auto div = gridCells[right].u - gridCells[center].u + gridCells[top].v - gridCells[center].v;

					auto k           = 1.0;
					auto compression = cell.particleDensity - particleRestDensity;
					if (compression > 0.0)
						div = div - k * compression;

					auto pValue = -div / s;
					pValue *= overRelaxation;

					gridCells[center].u -= sx0 * pValue;
					gridCells[right].u += sx1 * pValue;
					gridCells[center].v -= sy0 * pValue;
					gridCells[top].v += sy1 * pValue;
				}
			}
		};
		for_each(std::execution::seq, gridCells.begin(), gridCells.end(), parseIncompressibility);
	}
}