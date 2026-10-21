module;

#include <algorithm>
#include <execution>
#include <vector>

export module ParticleIntegration;

import BaseStructures;
import Constants;

// Positions are stored in pixels, so express gravity in pixels/s^2 by
// converting the usual cell-space acceleration to the grid's pixel scale.
constexpr float gravity = -9.81f * constants::cellHeight;

export void integrateParticles(std::vector<Particle>& particleMap, const float dt)
{
	auto integration = [dt](Particle& particle) {
		particle.velY += gravity * dt;
		particle.posX += particle.velX * dt;
		particle.posY += particle.velY * dt;
	};
	std::for_each(particleMap.begin(), particleMap.end(), integration);
}
