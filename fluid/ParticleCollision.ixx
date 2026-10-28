module;

#include <algorithm>
#include <ranges>
#include "ParallelForEach.h"

export module ParticleCollision;

import BaseStructures;
import Constants;
import Obstacle;

// can be modified at any time, they are not exported!! haha!
constexpr float minDist  = ObstacleConstants::radius + constants::particleRadius;
constexpr float minDist2 = minDist * minDist;

constexpr float LeftWall   = constants::cellHeight + constants::particleRadius;
constexpr float RightWall  = (constants::fNumX - 1) * constants::cellHeight - constants::particleRadius;
constexpr float bottomWall = constants::cellHeight + constants::particleRadius;
constexpr float topWall    = (constants::fNumY - 1) * constants::cellHeight - constants::particleRadius;


export void handleParticleCollisions(ParticleStore& particleMap, const Obstacle obstacle)
{
	auto particleCollision = [&](const std::size_t id) {
		const float x = particleMap.posX[id];
		const float y = particleMap.posY[id];
		const float dx = x - obstacle.x;
		const float dy = y - obstacle.y;
		const float d2 = dx * dx + dy * dy;

		// Match the reference FLIP demo: particles within the combined radius
		// inherit the obstacle velocity. The reference leaves their positions
		// untouched here; particle separation handles overlap independently.
		if (d2 < minDist2) {
			particleMap.velX[id] = obstacle.velX;
			particleMap.velY[id] = obstacle.velY;
		}

		// wall collisions
		bool hitsLeft  = x < LeftWall;
		bool hitsRight = x > RightWall;
		particleMap.posX[id] = std::ranges::clamp(x, LeftWall, RightWall);
		particleMap.velX[id] = hitsLeft || hitsRight ? 0.0f : particleMap.velX[id];

		bool hitsTop    = y > topWall;
		bool hitsBottom = y < bottomWall;
		particleMap.posY[id] = std::ranges::clamp(y, bottomWall, topWall);
		particleMap.velY[id] = hitsTop || hitsBottom ? 0.0f : particleMap.velY[id];
	};
	fluidParallelForEachIndex(particleMap.size(), particleCollision);
}
