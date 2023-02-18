module;

#include <algorithm>
#include <execution>
#include <vector>

export module ParticleIntegration;

import BaseStructures;
import Constants;

export void integrateParticles(std::vector<Particle>& particleMap)
{
	auto integration = [](Particle& particle) {
		particle.velY += constants::integ; //dt * gravity
		particle.posX += particle.velX * constants::dt;
		particle.posY += particle.velY * constants::dt;
	};
	std::for_each(std::execution::seq, particleMap.begin(), particleMap.end(), integration);
}