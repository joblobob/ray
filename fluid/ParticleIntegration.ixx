module;

#include <algorithm>
#include <execution>

export module ParticleIntegration;

import BaseStructures;
import Constants;

constexpr float gravity = -9.81f;
constexpr float integ   = constants::dt * gravity;

export void integrateParticles(std::vector<Particle>& particleMap)
{
	auto integration = [](Particle& particle) {
		particle.velY += particle.id < 3000 ? integ * 3.50f : particle.id < 15000 ? integ * 2.66f : integ; //dt * gravity * particle property magic :P
		particle.posX += particle.velX * constants::dt;
		particle.posY += particle.velY * constants::dt;
	};
	std::for_each(std::execution::seq, particleMap.begin(), particleMap.end(), integration);
}