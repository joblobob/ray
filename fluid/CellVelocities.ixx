module;

#include <algorithm>
#include <execution>

export module CellVelocities;

import BaseStructures;
import CellCalculations;


void restoreSolidCells(Cell& cell, const std::vector<Cell>& gridCells, const int fNumY)
{
	if (cell.du > 0.0)
		cell.u /= cell.du;

	if (cell.dv > 0.0)
		cell.v /= cell.dv;

	auto solid = cell.cellType == constants::CellType::Solid;
	if (solid || (cell.cellNumX > 0 && gridCells[(cell.cellNumX - 1) * fNumY + cell.cellNumY].cellType == constants::CellType::Solid))
		cell.u = cell.prevU;
	if (solid || (cell.cellNumY > 0 && gridCells[cell.cellNumX * fNumY + cell.cellNumY - 1].cellType == constants::CellType::Solid))
		cell.v = cell.prevV;
}

void setVelComponent(double& f, double& d, double pv, double delta)
{
	f += pv * delta;
	d += delta;
}

void parseVelocities(std::vector<Cell>& gridCells, Particle& particle, const double h, const int fNumX, const int fNumY, const double fInvSpacing)
{
	double h2 = 0.5 * h;
	int n     = fNumY;

	double dx = 0.0;
	double dy = h2;

	const double x = std::clamp(particle.posX, h, (fNumX - 1) * h);
	const double y = std::clamp(particle.posY, h, (fNumY - 1) * h);

	int x0    = std::min((int)floor((x - dx) * fInvSpacing), fNumX - 2);
	double tx = ((x - dx) - x0 * h) * fInvSpacing;
	int x1    = std::min(x0 + 1, fNumX - 2);

	int y0    = std::min((int)floor((y - dy) * fInvSpacing), fNumY - 2);
	double ty = ((y - dy) - y0 * h) * fInvSpacing;
	int y1    = std::min(y0 + 1, fNumY - 2);

	double sx = 1.0 - tx;
	double sy = 1.0 - ty;

	double d0 = sx * sy;
	double d1 = tx * sy;
	double d2 = tx * ty;
	double d3 = sx * ty;

	int nr0 = x0 * n + y0;
	int nr1 = x1 * n + y0;
	int nr2 = x1 * n + y1;
	int nr3 = x0 * n + y1;

	Cell& cell0 = gridCells[nr0];
	Cell& cell1 = gridCells[nr1];
	Cell& cell2 = gridCells[nr2];
	Cell& cell3 = gridCells[nr3];

	//specific

	setVelComponent(cell0.u, cell0.du, particle.velX, d0);
	setVelComponent(cell1.u, cell1.du, particle.velX, d1);
	setVelComponent(cell2.u, cell2.du, particle.velX, d2);
	setVelComponent(cell3.u, cell3.du, particle.velX, d3);


	dx = h2;
	dy = 0.0;

	x0 = std::min((int)floor((x - dx) * fInvSpacing), fNumX - 2);
	tx = ((x - dx) - x0 * h) * fInvSpacing;
	x1 = std::min(x0 + 1, fNumX - 2);

	y0 = std::min((int)floor((y - dy) * fInvSpacing), fNumY - 2);
	ty = ((y - dy) - y0 * h) * fInvSpacing;
	y1 = std::min(y0 + 1, fNumY - 2);

	sx = 1.0 - tx;
	sy = 1.0 - ty;

	d0 = sx * sy;
	d1 = tx * sy;
	d2 = tx * ty;
	d3 = sx * ty;

	nr0 = x0 * n + y0;
	nr1 = x1 * n + y0;
	nr2 = x1 * n + y1;
	nr3 = x0 * n + y1;

	Cell& cell0v = gridCells[nr0];
	Cell& cell1v = gridCells[nr1];
	Cell& cell2v = gridCells[nr2];
	Cell& cell3v = gridCells[nr3];

	setVelComponent(cell0v.v, cell0v.dv, particle.velY, d0);
	setVelComponent(cell1v.v, cell1v.dv, particle.velY, d1);
	setVelComponent(cell2v.v, cell2v.dv, particle.velY, d2);
	setVelComponent(cell3v.v, cell3v.dv, particle.velY, d3);
}



export void transferVelocitiesToGrid(std::vector<Particle>& particleMap,
    std::vector<Cell>& gridCells,
    const double h,
    const int fNumX,
    const int fNumY,
    const double fInvSpacing,
    const Border& fBorder)
{
	auto prepareCellType = [&](Cell& cell) {
		cell.prevU    = cell.u;
		cell.prevV    = cell.v;
		cell.du       = 0.0;
		cell.dv       = 0.0;
		cell.u        = 0.0;
		cell.v        = 0.0;
		cell.cellType = isVeryCloseToZero(cell.s) ? constants::CellType::Solid : constants::CellType::Air;
	};
	for_each(std::execution::seq, gridCells.begin(), gridCells.end(), prepareCellType);

	auto setCellType = [&](const Particle& p) {
		int cellNr     = cellNumber(p.posX, p.posY, fBorder, fInvSpacing);
		auto& cellType = gridCells[cellNr].cellType;
		if (cellType == constants::CellType::Air)
			cellType = constants::CellType::Fluid;
	};
	for_each(std::execution::seq, particleMap.begin(), particleMap.end(), setCellType);


	for_each(std::execution::seq, particleMap.begin(), particleMap.end(), [&](Particle& particle) {
		parseVelocities(gridCells, particle, h, fNumX, fNumY, fInvSpacing);
	});

	for_each(std::execution::seq, gridCells.begin(), gridCells.end(), [&](Cell& cell) { restoreSolidCells(cell, gridCells, fNumY); });
}
