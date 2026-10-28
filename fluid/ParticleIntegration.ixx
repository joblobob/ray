module;

#include <algorithm>
#include <vector>
#include "ParallelForEach.h"

export module ParticleIntegration;

import BaseStructures;
import Constants;

// Positions are stored in pixels. Assume the scene height represents 4 m,
// then convert physical gravity (m/s^2) to pixel acceleration (px/s^2).
// This is semi-implicit (symplectic) Euler: update velocity, then position.
constexpr float gravity = -9.81f * constants::maxheight / constants::worldHeightMeters;

export void integrateParticles(ParticleStore& particleMap, const float dt)
{
	auto integration = [&](const std::size_t id) {
		particleMap.velY[id] += gravity * dt;
		particleMap.posX[id] += particleMap.velX[id] * dt;
		particleMap.posY[id] += particleMap.velY[id] * dt;
	};
	fluidParallelForEachIndex(particleMap.size(), integration);
}
