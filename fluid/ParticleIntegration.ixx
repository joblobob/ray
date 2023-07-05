module;

#include <algorithm>
#include <ranges>
#include <vector>

export module ParticleIntegration;

import BaseStructures;
import Constants;

constexpr float gravity = -9.81f * 16.0f;
constexpr float integ   = constants::dt * gravity;

export void integrateParticles(std::vector<Particle>& particleMap)
{
	auto integration = [](Particle& particle) {
		particle.velY += particle.id < 3000  ? integ * 3.50f :
		                 particle.id < 10000 ? integ * 2.66f :
		                 particle.id < 15000 ? integ * 2.00f :
		                 particle.id < 17000 ? integ * 0.90f :
		                                       integ; //dt * gravity * particle property magic :P
		particle.posX += particle.velX * constants::dt;
		particle.posY += particle.velY * constants::dt;
		return particle;
	};
	particleMap | std::views::transform(integration) | std::ranges::to<std::vector>();
}