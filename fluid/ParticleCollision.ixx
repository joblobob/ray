module;

#include <algorithm>
#include <execution>
#include <ranges>

export module ParticleCollision;

import BaseStructures;
import Constants;
import Obstacle;

constexpr double minDist  = ObstacleConstants::radius + constants::particleRadius;
constexpr double minDist2 = minDist * minDist;

constexpr double minX = constants::cellHeight + constants::particleRadius;
constexpr double maxX = (constants::fNumX - 1) * constants::cellHeight - constants::particleRadius;
constexpr double minY = constants::cellHeight + constants::particleRadius;
constexpr double maxY = (constants::fNumY - 1) * constants::cellHeight - constants::particleRadius;

export void handleParticleCollisions(std::vector<Particle>& particleMap, const Obstacle obstacle)
{
	auto particleCollision = [&](Particle& particle) {
		double x = particle.posX;
		double y = particle.posY;

		const double dx = x - obstacle.x;
		const double dy = y - obstacle.y;
		const double d2 = dx * dx + dy * dy;

		// obstacle collision

		if (d2 < minDist2) {
			particle.velX = obstacle.velX;
			particle.velY = obstacle.velY;
		}
		// wall collisions

		if (x < minX) {
			x             = minX;
			particle.velX = 0.0;
		} else if (x > maxX) {
			x             = maxX;
			particle.velX = 0.0;
		}
		if (y < minY) {
			y             = minY;
			particle.velY = 0.0;
		} else if (y > maxY) {
			y             = maxY;
			particle.velY = 0.0;
		}
		particle.posX = x;
		particle.posY = y;
	};
	std::ranges::for_each(particleMap, particleCollision);
}