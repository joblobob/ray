module;

#include <algorithm>
#include <execution>

export module ParticleVelocities;

import BaseStructures;

//constants
constexpr double flipRatio = 0.9;

export void transferVelocitiesToParticles(std::vector<Particle>& particleMap, std::vector<Cell>& gridCells)
{
	const auto parseVelocitiesParticles = [&](Particle& particle) {
		double dx = 0.0;
		double dy = constants::halfCellHeight;

		const double x = std::clamp(particle.posX, constants::cellHeight, (constants::fNumX - 1) * constants::cellHeight);
		const double y = std::clamp(particle.posY, constants::cellHeight, (constants::fNumY - 1) * constants::cellHeight);

		int x0    = std::min(int_floor((x - dx) * constants::fInvSpacing), constants::fNumX - 2);
		double tx = ((x - dx) - x0 * constants::cellHeight) * constants::fInvSpacing;
		int x1    = std::min(x0 + 1, constants::fNumX - 2);

		int y0    = std::min(int_floor((y - dy) * constants::fInvSpacing), constants::fNumY - 2);
		double ty = ((y - dy) - y0 * constants::cellHeight) * constants::fInvSpacing;
		int y1    = std::min(y0 + 1, constants::fNumY - 2);

		double sx = 1.0 - tx;
		double sy = 1.0 - ty;

		double d0 = sx * sy;
		double d1 = tx * sy;
		double d2 = tx * ty;
		double d3 = sx * ty;

		int nr0 = x0 * constants::fNumY + y0;
		int nr1 = x1 * constants::fNumY + y0;
		int nr2 = x1 * constants::fNumY + y1;
		int nr3 = x0 * constants::fNumY + y1;

		Cell& cell0 = gridCells[nr0];
		Cell& cell1 = gridCells[nr1];
		Cell& cell2 = gridCells[nr2];
		Cell& cell3 = gridCells[nr3];

		// specific
		const Cell cellOffset0 = gridCells[nr0 - constants::fNumY];
		const Cell cellOffset1 = gridCells[nr1 - constants::fNumY];
		const Cell cellOffset2 = gridCells[nr2 - constants::fNumY];
		const Cell cellOffset3 = gridCells[nr3 - constants::fNumY];

		double valid0 = cell0.cellType != constants::CellType::Air || cellOffset0.cellType != constants::CellType::Air ? 1.0 : 0.0;
		double valid1 = cell1.cellType != constants::CellType::Air || cellOffset1.cellType != constants::CellType::Air ? 1.0 : 0.0;
		double valid2 = cell2.cellType != constants::CellType::Air || cellOffset2.cellType != constants::CellType::Air ? 1.0 : 0.0;
		double valid3 = cell3.cellType != constants::CellType::Air || cellOffset3.cellType != constants::CellType::Air ? 1.0 : 0.0;

		double v = particle.velX;
		double d = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

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
		dy = 0.0;

		x0 = std::min(int_floor((x - dx) * constants::fInvSpacing), constants::fNumX - 2);
		tx = ((x - dx) - x0 * constants::cellHeight) * constants::fInvSpacing;
		x1 = std::min(x0 + 1, constants::fNumX - 2);

		y0 = std::min(int_floor((y - dy) * constants::fInvSpacing), constants::fNumY - 2);
		ty = ((y - dy) - y0 * constants::cellHeight) * constants::fInvSpacing;
		y1 = std::min(y0 + 1, constants::fNumY - 2);

		sx = 1.0 - tx;
		sy = 1.0 - ty;

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
		const Cell cellOffset0v = gridCells[nr0 - 1];
		const Cell cellOffset1v = gridCells[nr1 - 1];
		const Cell cellOffset2v = gridCells[nr2 - 1];
		const Cell cellOffset3v = gridCells[nr3 - 1];

		valid0 = cell0.cellType != constants::CellType::Air || cellOffset0v.cellType != constants::CellType::Air ? 1.0 : 0.0;
		valid1 = cell1.cellType != constants::CellType::Air || cellOffset1v.cellType != constants::CellType::Air ? 1.0 : 0.0;
		valid2 = cell2.cellType != constants::CellType::Air || cellOffset2v.cellType != constants::CellType::Air ? 1.0 : 0.0;
		valid3 = cell3.cellType != constants::CellType::Air || cellOffset3v.cellType != constants::CellType::Air ? 1.0 : 0.0;

		v = particle.velY;
		d = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

		if (d > 0.0) {
			auto picV = (valid0 * d0 * cell0v.v + valid1 * d1 * cell1v.v + valid2 * d2 * cell2v.v + valid3 * d3 * cell3v.v) / d;
			auto corr = (valid0 * d0 * (cell0v.v - cell0v.prevV) + valid1 * d1 * (cell1v.v - cell1v.prevV) + valid2 * d2 * (cell2v.v - cell2v.prevV) +
			                valid3 * d3 * (cell3v.v - cell3v.prevV)) /
			            d;
			auto flipV = v + corr;

			auto newVelValue = (1.0 - flipRatio) * picV + flipRatio * flipV;
			particle.velY    = newVelValue;
		}
	};
	for_each(std::execution::par_unseq, particleMap.begin(), particleMap.end(), parseVelocitiesParticles);
}
