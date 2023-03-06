module;

#include <algorithm>
#include <execution>

export module PushParticles;

import BaseStructures;
import Constants;
import CellCalculations;

//constants
constexpr float pInvSpacing = 1.0f / (2.2f * constants::particleRadius);
constexpr int pNumX         = static_cast<int>((constants::maxwidth * pInvSpacing) + 1);
constexpr int pNumY         = static_cast<int>((constants::maxheight * pInvSpacing) + 1);
constexpr int pNumCells     = pNumX * pNumY;

constexpr int numParticleIters      = 2;
constexpr float colorDiffusionCoeff = 0.001f;
constexpr float minDist             = 2.0f * constants::particleRadius;
constexpr float minDist2            = minDist * minDist;

inline void calcColor(float& p1Color, float& p2Color)
{
	const float color = (p1Color + p2Color) * 0.5f;
	p1Color           = p1Color + (color - p1Color) * colorDiffusionCoeff;
	p2Color           = p2Color + (color - p2Color) * colorDiffusionCoeff;
}

static const Border pBorder { .maxX = pNumX, .maxY = pNumY };

export void pushParticlesApart(std::vector<Particle>& particleMap)
{
	std::vector<ParticleInCells> particleCells(pNumCells + 1, { 0, 0 }); // particlesincellinfo

	std::vector<int> cellNumbers(constants::maxParticles, 0);

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

	std::vector<int> cellParticleIds(constants::maxParticles, 0);

	// fill particles into cells
	auto fillParticleIntoCell = [&](const Particle& p) {
		particleCells[cellNumbers[p.id]].firstCellParticle--;
		cellParticleIds[particleCells[cellNumbers[p.id]].firstCellParticle] = p.id;
	};
	for_each(std::execution::seq, particleMap.begin(), particleMap.end(), fillParticleIntoCell);


	// push particles apart


	auto pushParticles = [&](Particle& particle) {
		const float px = particle.posX;
		const float py = particle.posY;

		const int pxi = static_cast<int>(px * pInvSpacing);
		const int pyi = static_cast<int>(py * pInvSpacing);
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
						const float qx = particleAtId.posX;
						const float qy = particleAtId.posY;

						float dx       = qx - px;
						float dy       = qy - py;
						const float d2 = dx * dx + dy * dy;
						if (!(d2 > minDist2 || d2 < 0.1)) {
							const float d = sqrt(d2);
							const float s = 0.5f * (minDist - d) / d;
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

	for (int iter = 0; iter < numParticleIters; iter++) {
		for_each(std::execution::par_unseq, particleMap.begin(), particleMap.end(), pushParticles);
	}
}