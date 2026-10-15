module;

#include <algorithm>
#include <execution>
#include <vector>

export module ParticleVelocities;

import BaseStructures;
import Constants;

//constants
constexpr float flipRatio = 0.9f;

export void transferVelocitiesToParticles(std::vector<Particle>& particleMap, std::vector<Cell>& gridCells)
{
	const auto parseVelocitiesParticles = [&](Particle& particle) {
		float dx = 0.0f;
		float dy = constants::halfCellHeight;

		const float x = std::clamp(particle.posX, constants::cellHeight, (constants::fNumX - 1) * constants::cellHeight);
		const float y = std::clamp(particle.posY, constants::cellHeight, (constants::fNumY - 1) * constants::cellHeight);

		int x0   = std::min(static_cast<int>((x - dx) * constants::fInvSpacing), constants::fNumX - 2);
		float tx = ((x - dx) - x0 * constants::cellHeight) * constants::fInvSpacing;
		int x1   = std::min(x0 + 1, constants::fNumX - 2);

		int y0   = std::min(static_cast<int>((y - dy) * constants::fInvSpacing), constants::fNumY - 2);
		float ty = ((y - dy) - y0 * constants::cellHeight) * constants::fInvSpacing;
		int y1   = std::min(y0 + 1, constants::fNumY - 2);

		float sx = 1.0f - tx;
		float sy = 1.0f - ty;

		float d0 = sx * sy;
		float d1 = tx * sy;
		float d2 = tx * ty;
		float d3 = sx * ty;

		int nr0 = x0 * constants::fNumY + y0;
		int nr1 = x1 * constants::fNumY + y0;
		int nr2 = x1 * constants::fNumY + y1;
		int nr3 = x0 * constants::fNumY + y1;

		Cell& cell0 = gridCells[nr0];
		Cell& cell1 = gridCells[nr1];
		Cell& cell2 = gridCells[nr2];
		Cell& cell3 = gridCells[nr3];

		// specific
		// A particle touching the left wall can sample x0 == 0. In that
		// case there is no cell one column to the left; use the boundary cell
		// itself, which is solid, instead of indexing before the grid.
		const auto leftCell = [&](const int cellNumber) -> const Cell& {
			return cellNumber >= constants::fNumY ? gridCells[cellNumber - constants::fNumY] : gridCells[cellNumber];
		};
		const Cell& cellOffset0 = leftCell(nr0);
		const Cell& cellOffset1 = leftCell(nr1);
		const Cell& cellOffset2 = leftCell(nr2);
		const Cell& cellOffset3 = leftCell(nr3);

		float valid0 = cell0.cellType != constants::CellType::Air || cellOffset0.cellType != constants::CellType::Air ? 1.0 : 0.0;
		float valid1 = cell1.cellType != constants::CellType::Air || cellOffset1.cellType != constants::CellType::Air ? 1.0 : 0.0;
		float valid2 = cell2.cellType != constants::CellType::Air || cellOffset2.cellType != constants::CellType::Air ? 1.0 : 0.0;
		float valid3 = cell3.cellType != constants::CellType::Air || cellOffset3.cellType != constants::CellType::Air ? 1.0 : 0.0;

		float v = particle.velX;
		float d = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

		if (d > 0.0) {
			auto picV = (valid0 * d0 * cell0.u + valid1 * d1 * cell1.u + valid2 * d2 * cell2.u + valid3 * d3 * cell3.u) / d;
			auto corr = (valid0 * d0 * (cell0.u - cell0.prevU) + valid1 * d1 * (cell1.u - cell1.prevU) + valid2 * d2 * (cell2.u - cell2.prevU) +
			                valid3 * d3 * (cell3.u - cell3.prevU)) /
			            d;
			auto flipV = v + corr;

			auto newVelValue = (1.0 - flipRatio) * picV + flipRatio * flipV;
			particle.velX    = newVelValue;
		}

		// 2eme passe
		dx = constants::halfCellHeight;
		dy = 0.0f;

		x0 = std::min(static_cast<int>((x - dx) * constants::fInvSpacing), constants::fNumX - 2);
		tx = ((x - dx) - x0 * constants::cellHeight) * constants::fInvSpacing;
		x1 = std::min(x0 + 1, constants::fNumX - 2);

		y0 = std::min(static_cast<int>((y - dy) * constants::fInvSpacing), constants::fNumY - 2);
		ty = ((y - dy) - y0 * constants::cellHeight) * constants::fInvSpacing;
		y1 = std::min(y0 + 1, constants::fNumY - 2);

		sx = 1.0f - tx;
		sy = 1.0f - ty;

		d0 = sx * sy;
		d1 = tx * sy;
		d2 = tx * ty;
		d3 = sx * ty;

		nr0 = x0 * constants::fNumY + y0;
		nr1 = x1 * constants::fNumY + y0;
		nr2 = x1 * constants::fNumY + y1;
		nr3 = x0 * constants::fNumY + y1;

		Cell& cell0v = gridCells[nr0];
		Cell& cell1v = gridCells[nr1];
		Cell& cell2v = gridCells[nr2];
		Cell& cell3v = gridCells[nr3];

		// specific
		// Likewise, y0 == 0 is valid at the bottom wall but has no row below
		// it. Treat the boundary cell as its own solid neighbour.
		const auto bottomCell = [&](const int cellNumber) -> const Cell& {
			return cellNumber % constants::fNumY > 0 ? gridCells[cellNumber - 1] : gridCells[cellNumber];
		};
		const Cell& cellOffset0v = bottomCell(nr0);
		const Cell& cellOffset1v = bottomCell(nr1);
		const Cell& cellOffset2v = bottomCell(nr2);
		const Cell& cellOffset3v = bottomCell(nr3);

		valid0 = cell0.cellType != constants::CellType::Air || cellOffset0v.cellType != constants::CellType::Air ? 1.0f : 0.0f;
		valid1 = cell1.cellType != constants::CellType::Air || cellOffset1v.cellType != constants::CellType::Air ? 1.0f : 0.0f;
		valid2 = cell2.cellType != constants::CellType::Air || cellOffset2v.cellType != constants::CellType::Air ? 1.0f : 0.0f;
		valid3 = cell3.cellType != constants::CellType::Air || cellOffset3v.cellType != constants::CellType::Air ? 1.0f : 0.0f;

		v = particle.velY;
		d = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

		if (d > 0.0f) {
			auto picV = (valid0 * d0 * cell0v.v + valid1 * d1 * cell1v.v + valid2 * d2 * cell2v.v + valid3 * d3 * cell3v.v) / d;
			auto corr = (valid0 * d0 * (cell0v.v - cell0v.prevV) + valid1 * d1 * (cell1v.v - cell1v.prevV) + valid2 * d2 * (cell2v.v - cell2v.prevV) +
			                valid3 * d3 * (cell3v.v - cell3v.prevV)) /
			            d;
			auto flipV = v + corr;

			auto newVelValue = (1.0f - flipRatio) * picV + flipRatio * flipV;
			particle.velY    = newVelValue;
		}
	};
	std::for_each(particleMap.begin(), particleMap.end(), parseVelocitiesParticles);
}
