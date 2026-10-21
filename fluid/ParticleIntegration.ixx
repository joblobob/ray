module;

#include <algorithm>
#include <vector>
#include "ParallelForEach.h"

export module ParticleIntegration;

import BaseStructures;
import Constants;

// Positions are stored in pixels. Assume the scene height represents 4 m,
// then convert physical gravity (m/s^2) to pixel acceleration (px/s^2).
constexpr float gravity = -9.81f * constants::maxheight / constants::worldHeightMeters;

export void integrateParticles(std::vector<Particle>& particleMap, const float dt)
{
	auto integration = [dt](Particle& particle) {
		particle.velY += gravity * dt;
		particle.posX += particle.velX * dt;
		particle.posY += particle.velY * dt;
	};
	fluidParallelForEach(particleMap.begin(), particleMap.end(), integration);
}
