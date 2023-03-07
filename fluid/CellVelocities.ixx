module;

#include <algorithm>
#include <execution>

export module CellVelocities;

import BaseStructures;
import CellCalculations;


void restoreCells(Cell& cell, const std::vector<Cell>& gridCells)
{
	if (cell.du > 0.0f)
		cell.u /= cell.du;

	if (cell.dv > 0.0f)
		cell.v /= cell.dv;

	if (isSolid(cell))
		restoreSolidCells(cell, gridCells);
}

void restoreSolidCells(Cell& cell, const std::vector<Cell>& gridCells)
{
	if (cell.cellNumX > 0 && gridCells[(cell.cellNumX - 1) * constants::fNumY + cell.cellNumY].cellType == constants::CellType::Solid)
		cell.u = cell.prevU;
	if (cell.cellNumY > 0 && gridCells[cell.cellNumX * constants::fNumY + cell.cellNumY - 1].cellType == constants::CellType::Solid)
		cell.v = cell.prevV;
}

void setVelComponent(float& f, float& d, const float pv, const float delta)
{
	f += pv * delta;
	d += delta;
}

void parseVelocitiesU(std::vector<Cell>& gridCells, const Particle& particle)
{
	const float x = std::clamp(particle.posX, constants::cellHeight, (constants::fNumX - 1) * constants::cellHeight);
	const float y = std::clamp(particle.posY, constants::cellHeight, (constants::fNumY - 1) * constants::cellHeight);

	const int x0   = std::min(static_cast<int>(x * constants::fInvSpacing), constants::fNumX - 2);
	const float tx = (x - x0 * constants::cellHeight) * constants::fInvSpacing;
	const int x1   = std::min(x0 + 1, constants::fNumX - 2);

	const int y0   = std::min(static_cast<int>((y - constants::halfCellHeight) * constants::fInvSpacing), constants::fNumY - 2);
	const float ty = ((y - constants::halfCellHeight) - y0 * constants::cellHeight) * constants::fInvSpacing;
	const int y1   = std::min(y0 + 1, constants::fNumY - 2);

	const float sx = 1.0f - tx;
	const float sy = 1.0f - ty;

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
	const float x = std::clamp(particle.posX, constants::cellHeight, (constants::fNumX - 1) * constants::cellHeight);
	const float y = std::clamp(particle.posY, constants::cellHeight, (constants::fNumY - 1) * constants::cellHeight);

	const int x0   = std::min(static_cast<int>((x - constants::halfCellHeight) * constants::fInvSpacing), constants::fNumX - 2);
	const float tx = ((x - constants::halfCellHeight) - x0 * constants::cellHeight) * constants::fInvSpacing;
	const int x1   = std::min(x0 + 1, constants::fNumX - 2);

	const int y0   = std::min(static_cast<int>((y - 0.0f) * constants::fInvSpacing), constants::fNumY - 2);
	const float ty = ((y - 0.0f) - y0 * constants::cellHeight) * constants::fInvSpacing;
	const int y1   = std::min(y0 + 1, constants::fNumY - 2);

	const float sx = 1.0 - tx;
	const float sy = 1.0 - ty;

	Cell& cell0v = gridCells[x0 * constants::fNumY + y0];
	Cell& cell1v = gridCells[x1 * constants::fNumY + y0];
	Cell& cell2v = gridCells[x1 * constants::fNumY + y1];
	Cell& cell3v = gridCells[x0 * constants::fNumY + y1];

	setVelComponent(cell0v.v, cell0v.dv, particle.velY, sx * sy);
	setVelComponent(cell1v.v, cell1v.dv, particle.velY, tx * sy);
	setVelComponent(cell2v.v, cell2v.dv, particle.velY, tx * ty);
	setVelComponent(cell3v.v, cell3v.dv, particle.velY, sx * ty);
}

void prepareCellType(Cell& cell)
{
	cell.prevU    = cell.u;
	cell.prevV    = cell.v;
	cell.du       = 0.0f;
	cell.dv       = 0.0f;
	cell.u        = 0.0f;
	cell.v        = 0.0f;
	cell.cellType = cell.s < 0.1f ? constants::CellType::Solid : constants::CellType::Air;
}

void setCellType(const Particle& p, std::vector<Cell>& gridCells, const Border& fBorder)
{
	int cellNr                    = cellNumber(p.posX, p.posY, fBorder, constants::fInvSpacing);
	constants::CellType& cellType = gridCells[cellNr].cellType;
	if (cellType == constants::CellType::Air)
		cellType = constants::CellType::Fluid;
}

export void transferVelocitiesToGrid(std::vector<Particle>& particleMap, std::vector<Cell>& gridCells, const Border& fBorder)
{
	for_each(std::execution::seq, gridCells.begin(), gridCells.end(), prepareCellType);

	for_each(std::execution::par_unseq, particleMap.cbegin(), particleMap.cend(), [&](const Particle& particle) {
		setCellType(particle, gridCells, fBorder);
		parseVelocitiesU(gridCells, particle);
		parseVelocitiesV(gridCells, particle);
	});

	for_each(std::execution::par_unseq, gridCells.begin(), gridCells.end(), restoreCells);
}
