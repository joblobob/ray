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

	const bool solid = cell.cellType == constants::CellType::Solid;
	if (solid || (cell.cellNumX > 0 && gridCells[(cell.cellNumX - 1) * fNumY + cell.cellNumY].cellType == constants::CellType::Solid))
		cell.u = cell.prevU;
	if (solid || (cell.cellNumY > 0 && gridCells[cell.cellNumX * fNumY + cell.cellNumY - 1].cellType == constants::CellType::Solid))
		cell.v = cell.prevV;
}

void setVelComponent(double& f, double& d, const double pv, const double delta)
{
	f += pv * delta;
	d += delta;
}

void parseVelocitiesU(std::vector<Cell>& gridCells, const Particle& particle)
{
	const double x = std::clamp(particle.posX, constants::cellHeight, (constants::fNumX - 1) * constants::cellHeight);
	const double y = std::clamp(particle.posY, constants::cellHeight, (constants::fNumY - 1) * constants::cellHeight);

	const int x0 = std::min(static_cast<int>(x * constants::fInvSpacing), constants::fNumX - 2);
	const double tx       = (x - x0 * constants::cellHeight) * constants::fInvSpacing;
	const int x1 = std::min(x0 + 1, constants::fNumX - 2);

	const int y0 = std::min(static_cast<int>((y - constants::halfCellHeight) * constants::fInvSpacing), constants::fNumY - 2);
	const double ty       = ((y - constants::halfCellHeight) - y0 * constants::cellHeight) * constants::fInvSpacing;
	const int y1 = std::min(y0 + 1, constants::fNumY - 2);

	const double sx = 1.0 - tx;
	const double sy = 1.0 - ty;

	Cell& cell0 = gridCells[x0 * constants::fNumY + y0];
	Cell& cell1 = gridCells[x1 * constants::fNumY + y0];
	Cell& cell2 = gridCells[x1 * constants::fNumY + y1];
	Cell& cell3 = gridCells[x0 * constants::fNumY + y1];

	//specific
	setVelComponent(cell0.u, cell0.du, particle.velX, sx * sy);
	setVelComponent(cell1.u, cell1.du, particle.velX, tx * sy);
	setVelComponent(cell2.u, cell2.du, particle.velX, tx * ty);
	setVelComponent(cell3.u, cell3.du, particle.velX, sx * ty);
}

void parseVelocitiesV(std::vector<Cell>& gridCells, const Particle& particle)
{
	const double x = std::clamp(particle.posX, constants::cellHeight, (constants::fNumX - 1) * constants::cellHeight);
	const double y = std::clamp(particle.posY, constants::cellHeight, (constants::fNumY - 1) * constants::cellHeight);

	const int x0 = std::min(static_cast<int>((x - constants::halfCellHeight) * constants::fInvSpacing), constants::fNumX - 2);
	const double tx       = ((x - constants::halfCellHeight) - x0 * constants::cellHeight) * constants::fInvSpacing;
	const int x1 = std::min(x0 + 1, constants::fNumX - 2);

	const int y0 = std::min(static_cast<int>((y - 0.0) * constants::fInvSpacing), constants::fNumY - 2);
	const double ty       = ((y - 0.0) - y0 * constants::cellHeight) * constants::fInvSpacing;
	const int y1 = std::min(y0 + 1, constants::fNumY - 2);

	const double sx = 1.0 - tx;
	const double sy = 1.0 - ty;

	Cell& cell0v = gridCells[x0 * constants::fNumY + y0];
	Cell& cell1v = gridCells[x1 * constants::fNumY + y0];
	Cell& cell2v = gridCells[x1 * constants::fNumY + y1];
	Cell& cell3v = gridCells[x0 * constants::fNumY + y1];

	setVelComponent(cell0v.v, cell0v.dv, particle.velY, sx * sy);
	setVelComponent(cell1v.v, cell1v.dv, particle.velY, tx * sy);
	setVelComponent(cell2v.v, cell2v.dv, particle.velY, tx * ty);
	setVelComponent(cell3v.v, cell3v.dv, particle.velY, sx * ty);
}



export void transferVelocitiesToGrid(std::vector<Particle>& particleMap, std::vector<Cell>& gridCells, const Border& fBorder)
{
	auto prepareCellType = [&](Cell& cell) {
		cell.prevU    = cell.u;
		cell.prevV    = cell.v;
		cell.du       = 0.0;
		cell.dv       = 0.0;
		cell.u        = 0.0;
		cell.v        = 0.0;
		cell.cellType = cell.s < 0.1 ? constants::CellType::Solid : constants::CellType::Air;
	};
	for_each(std::execution::par_unseq, gridCells.begin(), gridCells.end(), prepareCellType);

	auto setCellType = [&](const Particle& p) {
		int cellNr                    = cellNumber(p.posX, p.posY, fBorder, constants::fInvSpacing);
		constants::CellType& cellType = gridCells[cellNr].cellType;
		if (cellType == constants::CellType::Air)
			cellType = constants::CellType::Fluid;
	};
	for_each(std::execution::par_unseq, particleMap.cbegin(), particleMap.cend(), setCellType);


	for_each(std::execution::par_unseq, particleMap.cbegin(), particleMap.cend(), [&](const Particle& particle) {
		parseVelocitiesU(gridCells, particle);
		parseVelocitiesV(gridCells, particle);
	});

	for_each(std::execution::par_unseq, gridCells.begin(), gridCells.end(), [&](Cell& cell) { restoreSolidCells(cell, gridCells, constants::fNumY); });
}
