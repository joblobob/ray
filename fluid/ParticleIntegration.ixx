module;

#include <algorithm>
#include <execution>

export module ParticleIntegration;

import BaseStructures;
import Constants;

constexpr float gravity = -9.81f * 1.0f;
constexpr float integ   = constants::dt * gravity;

export void integrateParticles(std::vector<Particle>& particleMap)
{
	auto integration = [](Particle& particle) {
		particle.velY += integ; //dt * gravity
		particle.posX += particle.velX * constants::dt;
		particle.posY += particle.velY * constants::dt;
	};
	std::for_each(std::execution::seq, particleMap.begin(), particleMap.end(), integration);
}