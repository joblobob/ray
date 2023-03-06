module;

#include <algorithm>
#include <execution>
#include <ranges>

export module ParticleCollision;

import BaseStructures;
import Constants;
import Obstacle;

constexpr float minDist  = ObstacleConstants::radius + constants::particleRadius;
constexpr float minDist2 = minDist * minDist;

constexpr float minX = constants::cellHeight + constants::particleRadius;
constexpr float maxX = (constants::fNumX - 1) * constants::cellHeight - constants::particleRadius;
constexpr float minY = constants::cellHeight + constants::particleRadius;
constexpr float maxY = (constants::fNumY - 1) * constants::cellHeight - constants::particleRadius;

export void handleParticleCollisions(std::vector<Particle>& particleMap, const Obstacle obstacle)
{
	auto particleCollision = [&](Particle& particle) {
		float x = particle.posX;
		float y = particle.posY;

		const float dx = x - obstacle.x;
		const float dy = y - obstacle.y;
		const float d2 = dx * dx + dy * dy;

		// obstacle collision

		if (d2 < minDist2) {
			particle.velX = obstacle.velX;
			particle.velY = obstacle.velY;
		}
		// wall collisions

		if (x < minX) {
			x             = minX;
			particle.velX = 0.0f;
		} else if (x > maxX) {
			x             = maxX;
			particle.velX = 0.0f;
		}
		if (y < minY) {
			y             = minY;
			particle.velY = 0.0f;
		} else if (y > maxY) {
			y             = maxY;
			particle.velY = 0.0f;
		}
		particle.posX = x;
		particle.posY = y;
	};
	std::ranges::for_each(particleMap, particleCollision);
}