module;

#include <algorithm>
#include <ranges>
#include <cmath>
#include "ParallelForEach.h"

export module ParticleCollision;

import BaseStructures;
import Constants;
import Obstacle;

// can be modified at any time, they are not exported!! haha!
constexpr float minDist  = ObstacleConstants::radius + constants::particleRadius;
constexpr float minDist2 = minDist * minDist;
constexpr float obstacleVelocityTransfer = 0.75f;

constexpr float LeftWall   = constants::cellHeight + constants::particleRadius;
constexpr float RightWall  = (constants::fNumX - 1) * constants::cellHeight - constants::particleRadius;
constexpr float bottomWall = constants::cellHeight + constants::particleRadius;
constexpr float topWall    = (constants::fNumY - 1) * constants::cellHeight - constants::particleRadius;


export void handleParticleCollisions(std::vector<Particle>& particleMap, const Obstacle obstacle)
{
	auto particleCollision = [&](Particle& particle) {
		float x = particle.posX;
		float y = particle.posY;
		const float dx = x - obstacle.x;
		const float dy = y - obstacle.y;
		const float d2 = dx * dx + dy * dy;

		// Keep particle centers outside the obstacle and transfer only part of
		// its inward normal motion. Copying its full velocity to every overlap
		// creates a large, coherent kick and drains particles from the fluid.
		if (d2 < minDist2) {
			const float distance = std::sqrt(std::max(d2, 1.0e-6f));
			const float nx = d2 > 1.0e-6f ? dx / distance : 1.0f;
			const float ny = d2 > 1.0e-6f ? dy / distance : 0.0f;
			particle.posX = obstacle.x + nx * minDist;
			particle.posY = obstacle.y + ny * minDist;

			const float relativeNormalVelocity = (particle.velX - obstacle.velX) * nx + (particle.velY - obstacle.velY) * ny;
			if (relativeNormalVelocity < 0.0f) {
				const float impulse = -relativeNormalVelocity * obstacleVelocityTransfer;
				particle.velX += impulse * nx;
				particle.velY += impulse * ny;
			}
			x = particle.posX;
			y = particle.posY;
		}

		// wall collisions
		bool hitsLeft  = x < LeftWall;
		bool hitsRight = x > RightWall;
		particle.posX  = std::ranges::clamp(x, LeftWall, RightWall);
		particle.velX  = hitsLeft || hitsRight ? 0.0f : particle.velX;

		bool hitsTop    = y > topWall;
		bool hitsBottom = y < bottomWall;
		particle.posY   = std::ranges::clamp(y, bottomWall, topWall);
		particle.velY   = hitsTop || hitsBottom ? 0.0f : particle.velY;
	};
	fluidParallelForEach(particleMap.begin(), particleMap.end(), particleCollision);
}
