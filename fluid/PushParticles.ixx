module;

#include "ParallelForEach.h"
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <vector>

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
constexpr float halfColorDiffusion  = 0.5f * colorDiffusionCoeff;

export struct ParticleInCells {
	int numCellParticles = 0, firstCellParticle = 0;
};

export struct PushParticlesBuffers {
	std::vector<ParticleInCells> particleCells;
	std::vector<int> cellNumbers;
	std::vector<int> cellParticleIds;
	std::vector<std::array<float, 5> > particleCorrections;
};

static const Border pBorder { .maxX = pNumX, .maxY = pNumY };

export void pushParticlesApart(ParticleStore& particleMap, PushParticlesBuffers& buffers)
{
	buffers.particleCells.assign(pNumCells + 1, { 0, 0 });
	buffers.cellNumbers.resize(constants::maxParticles);
	buffers.cellParticleIds.resize(constants::maxParticles);

	// Phase 1: count particles per cell (parallel with atomics)
	const auto countParticlesInCell = [&](const int id) {
		const int cellnum       = cellNumber(particleMap.posX[id], particleMap.posY[id], pBorder, pInvSpacing);
		buffers.cellNumbers[id] = cellnum;
		std::atomic_ref<int> counter(buffers.particleCells[cellnum].numCellParticles);
		counter.fetch_add(1, std::memory_order_relaxed);
	};
	fluidParallelForEachIndex(particleMap.size(), countParticlesInCell);

	// Phase 2: exclusive prefix sum over grid cells (serial, fast - only ~4000 cells)
	int first = 0;
	for (auto& pIncell : buffers.particleCells) {
		const int count           = pIncell.numCellParticles;
		pIncell.firstCellParticle = first;
		first += count;
	}
	buffers.particleCells[pNumCells].firstCellParticle = first; // guard

	// Phase 3: fill particles into cells (parallel with atomics)
	const auto fillParticleIntoCell = [&](const int id) {
		const int cellnum = buffers.cellNumbers[id];
		std::atomic_ref<int> counter(buffers.particleCells[cellnum].firstCellParticle);
		const int idx                = counter.fetch_sub(1, std::memory_order_relaxed) - 1;
		buffers.cellParticleIds[idx] = id;
	};
	fluidParallelForEachIndex(particleMap.size(), fillParticleIntoCell);

	// Phase 4 & 5: gather and apply corrections (already parallel)
	buffers.particleCorrections.resize(particleMap.size());
	const auto gatherParticleCorrection = [&](const std::size_t particleIndex) {
		const float px    = particleMap.posX[particleIndex];
		const float py    = particleMap.posY[particleIndex];
		float correctionX = 0.0f;
		float correctionY = 0.0f;
		float correctionR = 0.0f;
		float correctionG = 0.0f;
		float correctionB = 0.0f;

		const int pxi = static_cast<int>(px * pInvSpacing);
		const int pyi = static_cast<int>(py * pInvSpacing);
		const int x0  = std::max(pxi - 1, 0);
		const int y0  = std::max(pyi - 1, 0);
		const int x1  = std::min(pxi + 1, pNumX - 1);
		const int y1  = std::min(pyi + 1, pNumY - 1);

		for (int xi = x0; xi <= x1; xi++) {
			for (int yi = y0; yi <= y1; yi++) {
				const int cellNr = xi * pNumY + yi;
				const int first  = buffers.particleCells[cellNr].firstCellParticle;
				const int last   = buffers.particleCells[cellNr + 1].firstCellParticle;
				for (int j = first; j < last; j++) {
					const int neighborId = buffers.cellParticleIds[j];
					if (neighborId == particleIndex)
						continue;

					const float dx              = particleMap.posX[neighborId] - px;
					const float dy              = particleMap.posY[neighborId] - py;
					const float distanceSquared = dx * dx + dy * dy;
					if (distanceSquared > minDist2 || distanceSquared < 0.1f)
						continue;

					const float distance = std::sqrt(distanceSquared);
					const float scale    = 0.5f * (minDist - distance) / distance;
					correctionX -= dx * scale;
					correctionY -= dy * scale;
					correctionR += (particleMap.colorR[neighborId] - particleMap.colorR[particleIndex]) * halfColorDiffusion;
					correctionG += (particleMap.colorG[neighborId] - particleMap.colorG[particleIndex]) * halfColorDiffusion;
					correctionB += (particleMap.colorB[neighborId] - particleMap.colorB[particleIndex]) * halfColorDiffusion;
				}
			}
		}
		buffers.particleCorrections[particleIndex] = { correctionX, correctionY, correctionR, correctionG, correctionB };
	};
	const auto applyParticleCorrection = [&](const std::size_t id) {
		const auto& correction = buffers.particleCorrections[id];
		particleMap.posX[id] += correction[0];
		particleMap.posY[id] += correction[1];
		particleMap.colorR[id] += correction[2];
		particleMap.colorG[id] += correction[3];
		particleMap.colorB[id] += correction[4];
	};

	for (int iter = 0; iter < numParticleIters; iter++) {
		fluidParallelForEachIndex(particleMap.size(), gatherParticleCorrection);
		fluidParallelForEachIndex(particleMap.size(), applyParticleCorrection);
	}
}