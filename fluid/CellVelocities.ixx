module;

#include "ParallelForEach.h"
#include <algorithm>
#include <atomic>
#include <vector>

export module CellVelocities;

import BaseStructures;
import CellCalculations;
import Constants;


void restoreSolidCells(Cell& cell, const std::vector<Cell>& gridCells)
{
	if (cell.cellNumX > 0 && gridCells[(cell.cellNumX - 1) * constants::fNumY + cell.cellNumY].cellType == constants::CellType::Solid)
		cell.u = cell.prevU;
	if (cell.cellNumY > 0 && gridCells[cell.cellNumX * constants::fNumY + cell.cellNumY - 1].cellType == constants::CellType::Solid)
		cell.v = cell.prevV;
}

void restoreCells(Cell& cell, const std::vector<Cell>& gridCells)
{
	if (cell.du > 0.0f)
		cell.u /= cell.du;

	if (cell.dv > 0.0f)
		cell.v /= cell.dv;

	if (isSolid(cell))
		restoreSolidCells(cell, gridCells);
}

void setVelComponent(float& f, float& d, const float pv, const float delta, const float particleMass)
{
	const float massWeight = delta * particleMass;
	f += pv * massWeight;
	d += massWeight;
}

void parseVelocitiesU(std::vector<Cell>& gridCells, const ParticleStore& particleMap, const int id)
{
	const float particleMass = constants::particleDensityRatioForId(id);
	const float x            = std::clamp(particleMap.posX[id], constants::cellHeight, (constants::fNumX - 1) * constants::cellHeight);
	const float y            = std::clamp(particleMap.posY[id], constants::cellHeight, (constants::fNumY - 1) * constants::cellHeight);

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
	setVelComponent(cell0.u, cell0.du, particleMap.velX[id], sx * sy, particleMass);
	setVelComponent(cell1.u, cell1.du, particleMap.velX[id], tx * sy, particleMass);
	setVelComponent(cell2.u, cell2.du, particleMap.velX[id], tx * ty, particleMass);
	setVelComponent(cell3.u, cell3.du, particleMap.velX[id], sx * ty, particleMass);
}

void parseVelocitiesV(std::vector<Cell>& gridCells, const ParticleStore& particleMap, const int id)
{
	const float particleMass = constants::particleDensityRatioForId(id);
	const float x            = std::clamp(particleMap.posX[id], constants::cellHeight, (constants::fNumX - 1) * constants::cellHeight);
	const float y            = std::clamp(particleMap.posY[id], constants::cellHeight, (constants::fNumY - 1) * constants::cellHeight);

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

	setVelComponent(cell0v.v, cell0v.dv, particleMap.velY[id], sx * sy, particleMass);
	setVelComponent(cell1v.v, cell1v.dv, particleMap.velY[id], tx * sy, particleMass);
	setVelComponent(cell2v.v, cell2v.dv, particleMap.velY[id], tx * ty, particleMass);
	setVelComponent(cell3v.v, cell3v.dv, particleMap.velY[id], sx * ty, particleMass);
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

void setCellType(const ParticleStore& particleMap, const int id, std::vector<Cell>& gridCells, const Border& fBorder)
{
	int cellNr                    = cellNumber(particleMap.posX[id], particleMap.posY[id], fBorder, constants::fInvSpacing);
	constants::CellType& cellType = gridCells[cellNr].cellType;
	if (cellType == constants::CellType::Air)
		cellType = constants::CellType::Fluid;
}

export void transferVelocitiesToGrid(const ParticleStore& particleMap, std::vector<Cell>& gridCells, const Border& fBorder)
{
	fluidParallelForEach(gridCells.begin(), gridCells.end(), prepareCellType);

	const auto processParticle = [&](const std::size_t id) {
		// setCellType
		const int cellNr = cellNumber(particleMap.posX[id], particleMap.posY[id], fBorder, constants::fInvSpacing);
		std::atomic_ref<constants::CellType> cellTypeRef(gridCells[cellNr].cellType);
		constants::CellType expected = constants::CellType::Air;
		cellTypeRef.compare_exchange_strong(expected, constants::CellType::Fluid, std::memory_order_relaxed);

		// parseVelocitiesU
		const float particleMass = constants::particleDensityRatioForId(static_cast<int>(id));
		const float x            = std::clamp(particleMap.posX[id], constants::cellHeight, (constants::fNumX - 1) * constants::cellHeight);
		const float y            = std::clamp(particleMap.posY[id], constants::cellHeight, (constants::fNumY - 1) * constants::cellHeight);

		const int x0   = std::min(static_cast<int>(x * constants::fInvSpacing), constants::fNumX - 2);
		const float tx = (x - x0 * constants::cellHeight) * constants::fInvSpacing;
		const int x1   = std::min(x0 + 1, constants::fNumX - 2);

		const int y0   = std::min(static_cast<int>((y - constants::halfCellHeight) * constants::fInvSpacing), constants::fNumY - 2);
		const float ty = ((y - constants::halfCellHeight) - y0 * constants::cellHeight) * constants::fInvSpacing;
		const int y1   = std::min(y0 + 1, constants::fNumY - 2);

		const float sx = 1.0f - tx;
		const float sy = 1.0f - ty;

		// U component
		Cell& cell0 = gridCells[x0 * constants::fNumY + y0];
		Cell& cell1 = gridCells[x1 * constants::fNumY + y0];
		Cell& cell2 = gridCells[x1 * constants::fNumY + y1];
		Cell& cell3 = gridCells[x0 * constants::fNumY + y1];

		const float pv = particleMap.velX[id];
		std::atomic_ref<float> du0(cell0.du);
		std::atomic_ref<float> u0(cell0.u);
		std::atomic_ref<float> du1(cell1.du);
		std::atomic_ref<float> u1(cell1.u);
		std::atomic_ref<float> du2(cell2.du);
		std::atomic_ref<float> u2(cell2.u);
		std::atomic_ref<float> du3(cell3.du);
		std::atomic_ref<float> u3(cell3.u);

		const float w0 = sx * sy;
		const float w1 = tx * sy;
		const float w2 = tx * ty;
		const float w3 = sx * ty;

		const float massWeight0 = w0 * particleMass;
		const float massWeight1 = w1 * particleMass;
		const float massWeight2 = w2 * particleMass;
		const float massWeight3 = w3 * particleMass;

		du0.fetch_add(massWeight0, std::memory_order_relaxed);
		u0.fetch_add(pv * massWeight0, std::memory_order_relaxed);

		du1.fetch_add(massWeight1, std::memory_order_relaxed);
		u1.fetch_add(pv * massWeight1, std::memory_order_relaxed);

		du2.fetch_add(massWeight2, std::memory_order_relaxed);
		u2.fetch_add(pv * massWeight2, std::memory_order_relaxed);

		du3.fetch_add(massWeight3, std::memory_order_relaxed);
		u3.fetch_add(pv * massWeight3, std::memory_order_relaxed);

		// V component
		const int x0v   = std::min(static_cast<int>((x - constants::halfCellHeight) * constants::fInvSpacing), constants::fNumX - 2);
		const float txv = ((x - constants::halfCellHeight) - x0v * constants::cellHeight) * constants::fInvSpacing;
		const int x1v   = std::min(x0v + 1, constants::fNumX - 2);

		const int y0v   = std::min(static_cast<int>((y - 0.0f) * constants::fInvSpacing), constants::fNumY - 2);
		const float tyv = ((y - 0.0f) - y0v * constants::cellHeight) * constants::fInvSpacing;
		const int y1v   = std::min(y0v + 1, constants::fNumY - 2);

		const float sxv = 1.0f - txv;
		const float syv = 1.0f - tyv;

		Cell& cell0v = gridCells[x0v * constants::fNumY + y0v];
		Cell& cell1v = gridCells[x1v * constants::fNumY + y0v];
		Cell& cell2v = gridCells[x1v * constants::fNumY + y1v];
		Cell& cell3v = gridCells[x0v * constants::fNumY + y1v];

		const float pvY = particleMap.velY[id];

		std::atomic_ref<float> dv0(cell0v.dv);
		std::atomic_ref<float> v0(cell0v.v);
		std::atomic_ref<float> dv1(cell1v.dv);
		std::atomic_ref<float> v1(cell1v.v);
		std::atomic_ref<float> dv2(cell2v.dv);
		std::atomic_ref<float> v2(cell2v.v);
		std::atomic_ref<float> dv3(cell3v.dv);
		std::atomic_ref<float> v3(cell3v.v);

		const float w0v = sxv * syv;
		const float w1v = txv * syv;
		const float w2v = txv * tyv;
		const float w3v = sxv * tyv;

		const float massWeight0v = w0v * particleMass;
		const float massWeight1v = w1v * particleMass;
		const float massWeight2v = w2v * particleMass;
		const float massWeight3v = w3v * particleMass;

		dv0.fetch_add(massWeight0v, std::memory_order_relaxed);
		v0.fetch_add(pvY * massWeight0v, std::memory_order_relaxed);

		dv1.fetch_add(massWeight1v, std::memory_order_relaxed);
		v1.fetch_add(pvY * massWeight1v, std::memory_order_relaxed);

		dv2.fetch_add(massWeight2v, std::memory_order_relaxed);
		v2.fetch_add(pvY * massWeight2v, std::memory_order_relaxed);

		dv3.fetch_add(massWeight3v, std::memory_order_relaxed);
		v3.fetch_add(pvY * massWeight3v, std::memory_order_relaxed);
	};
	fluidParallelForEachIndex(particleMap.size(), processParticle);

	fluidParallelForEach(gridCells.begin(), gridCells.end(), [&](Cell& cell) { restoreCells(cell, gridCells); });
}