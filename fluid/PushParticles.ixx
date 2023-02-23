module;

#include <algorithm>
#include <execution>

export module PushParticles;

import BaseStructures;
import Constants;
import CellCalculations;


inline void calcColor(double& p1Color, double& p2Color)
{
	const auto color = (p1Color + p2Color) * 0.5;
	p1Color          = p1Color + (color - p1Color) * constants::colorDiffusionCoeff;
	p2Color          = p2Color + (color - p2Color) * constants::colorDiffusionCoeff;
}

static const Border pBorder { .maxX = constants::pNumX, .maxY = constants::pNumY };

export void pushParticlesApart(std::vector<Particle>& particleMap,
    const int numIters,
    const int pNumCells,
    const int pNumX,
    const int pNumY,
    const int maxParticles,
    const double pInvSpacing,
    const double particleRadius)
{
	std::vector<ParticleInCells> particleCells(pNumCells + 1, { 0, 0 }); // particlesincellinfo

	std::vector<int> cellNumbers(maxParticles, 0);

	// count particles per cell
	const auto countParticlesInCell = [&](const Particle& particle) {
		const int cellnum        = cellNumber(particle.posX, particle.posY, pBorder, pInvSpacing);
		cellNumbers[particle.id] = cellnum;
		particleCells[cellnum].numCellParticles++;
	};
	for_each(std::execution::seq, particleMap.begin(), particleMap.end(), countParticlesInCell);

	// partial sums
	int first = 0;

	//count first cell
	const auto countFirstCell = [&](ParticleInCells& pIncell) {
		first += pIncell.numCellParticles;
		pIncell.firstCellParticle = first;
	};
	for_each(std::execution::seq, particleCells.begin(), particleCells.end(), countFirstCell);


	particleCells[pNumCells].firstCellParticle = first; // guard

	std::vector<int> cellParticleIds(maxParticles, 0);

	// fill particles into cells
	auto fillParticleIntoCell = [&](const Particle& p) {
		particleCells[cellNumbers[p.id]].firstCellParticle--;
		cellParticleIds[particleCells[cellNumbers[p.id]].firstCellParticle] = p.id;
	};
	for_each(std::execution::seq, particleMap.begin(), particleMap.end(), fillParticleIntoCell);


	// push particles apart
	const auto minDist  = 2.0 * particleRadius;
	const auto minDist2 = minDist * minDist;

	auto pushParticles = [&](auto& particle) {
		const auto px = particle.posX;
		const auto py = particle.posY;

		const int pxi = floor(px * pInvSpacing);
		const int pyi = floor(py * pInvSpacing);
		const int x0  = std::max(pxi - 1, 0);
		const int y0  = std::max(pyi - 1, 0);
		const int x1  = std::min(pxi + 1, pNumX - 1);
		const int y1  = std::min(pyi + 1, pNumY - 1);

		for (int xi = x0; xi <= x1; xi++) {
			for (int yi = y0; yi <= y1; yi++) {
				const int cellNr = xi * pNumY + yi;
				const int first  = particleCells[cellNr].firstCellParticle;
				const int last   = particleCells[cellNr + 1].firstCellParticle;
				for (int j = first; j < last; j++) {
					const int id = cellParticleIds[j];
					if (id != particle.id) {
						Particle& particleAtId { particleMap[id] };
						const double qx = particleAtId.posX;
						const double qy = particleAtId.posY;

						double dx     = qx - px;
						double dy     = qy - py;
						const auto d2 = dx * dx + dy * dy;
						if (!(d2 > minDist2 || isVeryCloseToZero(d2))) {
							const auto d = sqrt(d2);
							const auto s = 0.5 * (minDist - d) / d;
							dx *= s;
							dy *= s;
							particle.posX -= dx;
							particle.posY -= dy;
							particleAtId.posX += dx;
							particleAtId.posY += dy;

							// diffuse colors

							calcColor(particle.colorR, particleAtId.colorR);
							calcColor(particle.colorG, particleAtId.colorG);
							calcColor(particle.colorB, particleAtId.colorB);
						}
					}
				}
			}
		}
	};

	for (int iter = 0; iter < numIters; iter++) {
		for_each(std::execution::par_unseq, particleMap.begin(), particleMap.end(), pushParticles);
	}
}