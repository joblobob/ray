module;

#include <algorithm>
#include <execution>

export module ParticleVelocities;

import BaseStructures;


export void transferVelocitiesToParticles(std::vector<Particle>& particleMap,
    std::vector<Cell>& gridCells,
    const double h,
    const int fNumX,
    const int fNumY,
    const double fInvSpacing)
{
	const auto parseVelocitiesParticles = [&](Particle& particle) {
		const double h2 = 0.5 * h;
		const int n     = fNumY;

		double dx = 0.0;
		double dy = h2;

		const double x = std::clamp(particle.posX, h, (fNumX - 1) * h);
		const double y = std::clamp(particle.posY, h, (fNumY - 1) * h);

		int x0    = std::min(int_floor((x - dx) * fInvSpacing), fNumX - 2);
		double tx = ((x - dx) - x0 * h) * fInvSpacing;
		int x1    = std::min(x0 + 1, fNumX - 2);

		int y0    = std::min(int_floor((y - dy) * fInvSpacing), fNumY - 2);
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

		// specific

		int offset             = n;
		const Cell cellOffset0 = gridCells[nr0 - offset];
		const Cell cellOffset1 = gridCells[nr1 - offset];
		const Cell cellOffset2 = gridCells[nr2 - offset];
		const Cell cellOffset3 = gridCells[nr3 - offset];

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

			auto newVelValue = (1.0 - constants::flipRatio) * picV + constants::flipRatio * flipV;
			particle.velX    = newVelValue;
		}

		// 2eme passe
		dx = h2;
		dy = 0.0;

		x0 = std::min(int_floor((x - dx) * fInvSpacing), fNumX - 2);
		tx = ((x - dx) - x0 * h) * fInvSpacing;
		x1 = std::min(x0 + 1, fNumX - 2);

		y0 = std::min(int_floor((y - dy) * fInvSpacing), fNumY - 2);
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

		// specific


		offset                  = 1;
		const Cell cellOffset0v = gridCells[nr0 - offset];
		const Cell cellOffset1v = gridCells[nr1 - offset];
		const Cell cellOffset2v = gridCells[nr2 - offset];
		const Cell cellOffset3v = gridCells[nr3 - offset];

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

			auto newVelValue = (1.0 - constants::flipRatio) * picV + constants::flipRatio * flipV;
			particle.velY    = newVelValue;
		}
	};
	for_each(std::execution::par_unseq, particleMap.begin(), particleMap.end(), parseVelocitiesParticles);
}
