module;

#include <algorithm>
#include <execution>

export module ParticleCollision;

import BaseStructures;


export void handleParticleCollisions(std::vector<Particle>& particleMap,
    const int fNumX,
    const int fNumY,
    const double fInvSpacing,
    const double obstacleX,
    const double obstacleY,
    const double obstacleVelX,
    const double obstacleVelY,
    const double obstacleRadius,
    const double particleRadius)
{
	const auto h        = 1.0 / fInvSpacing;
	const auto r        = particleRadius;
	const auto minDist  = obstacleRadius + r;
	const auto minDist2 = minDist * minDist;

	const auto minX = h + r;
	const auto maxX = (fNumX - 1) * h - r;
	const auto minY = h + r;
	const auto maxY = (fNumY - 1) * h - r;

	auto particleCollision = [&](auto& particle) {
		auto x = particle.posX;
		auto y = particle.posY;

		const auto dx = x - obstacleX;
		const auto dy = y - obstacleY;
		const auto d2 = dx * dx + dy * dy;

		// obstacle collision

		if (d2 < minDist2) {
			particle.velX = obstacleVelX;
			particle.velY = obstacleVelY;
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
	for_each(std::execution::seq, particleMap.begin(), particleMap.end(), particleCollision);
}