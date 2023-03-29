module;

#include <algorithm>
#include <execution>
#include <ranges>

export module ParticleCollision;

import BaseStructures;
import Constants;
import Obstacle;

// can be modified at any time, they are not exported!! haha!


constexpr float LeftWall   = constants::cellHeight;
constexpr float RightWall  = (constants::fNumX - 1) * constants::cellHeight;
constexpr float bottomWall = constants::cellHeight;
constexpr float topWall    = (constants::fNumY - 1) * constants::cellHeight;


export void handleParticleCollisions(std::vector<Particle>& particleMap, const Obstacle obstacle)
{
	auto particleCollision = [&](Particle& particle) {
		const float minDist  = ObstacleConstants::radius + particle.radius;
		const float minDist2 = minDist * minDist;
		const float x        = particle.posX;
		const float y        = particle.posY;
		const float dx       = x - obstacle.x;
		const float dy       = y - obstacle.y;
		const float d2       = dx * dx + dy * dy;

		// obstacle collision
		bool hitsObstacle = d2 < minDist2;
		particle.velX     = hitsObstacle ? obstacle.velX : particle.velX;
		particle.velY     = hitsObstacle ? obstacle.velY : particle.velY;


		// wall collisions
		bool hitsLeft  = x < (LeftWall + particle.radius);
		bool hitsRight = x > (RightWall - particle.radius);
		particle.posX  = std::ranges::clamp(x, LeftWall, RightWall);
		particle.velX  = hitsLeft || hitsRight ? 0.0f : particle.velX;

		bool hitsTop    = y > (topWall - particle.radius);
		bool hitsBottom = y < (bottomWall + particle.radius);
		particle.posY   = std::ranges::clamp(y, bottomWall, topWall);
		particle.velY   = hitsTop || hitsBottom ? 0.0f : particle.velY;
	};
	std::for_each(std::execution::par_unseq, particleMap.begin(), particleMap.end(), particleCollision);
}