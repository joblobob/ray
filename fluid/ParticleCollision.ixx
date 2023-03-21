module;

#include <algorithm>
#include <execution>
#include <ranges>

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


export void handleParticleCollisions(std::vector<Particle>& particleMap, const Obstacle obstacle)
{
	auto particleCollision = [&](Particle& particle) {
		const float x  = particle.posX;
		const float y  = particle.posY;
		const float dx = x - obstacle.x;
		const float dy = y - obstacle.y;
		const float d2 = dx * dx + dy * dy;

		// obstacle collision
		bool hitsObstacle = [&]() {
			return d2 < minDist2;
		}();
		particle.velX = hitsObstacle ? obstacle.velX : particle.velX;
		particle.velY = hitsObstacle ? obstacle.velY : particle.velY;


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
	std::for_each(std::execution::par_unseq, particleMap.begin(), particleMap.end(), particleCollision);
}