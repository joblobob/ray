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

constexpr float LeftWall   = constants::cellHeight + constants::particleRadius;
constexpr float RightWall  = (constants::fNumX - 1) * constants::cellHeight - constants::particleRadius;
constexpr float topWall    = constants::cellHeight + constants::particleRadius;
constexpr float bottomWall = (constants::fNumY - 1) * constants::cellHeight - constants::particleRadius;


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
		bool hitsLeft = [&]() {
			return x < LeftWall;
		}();
		bool hitsRight = [&]() {
			return x > RightWall;
		}();
		particle.posX = hitsLeft ? LeftWall : hitsRight ? RightWall : x;
		particle.velX = hitsLeft ? 0.0f : hitsRight ? 0.0f : particle.velX;

		// yeah.. upside down
		bool hitsUp = [&]() {
			return y < topWall;
		}();
		bool hitsDown = [&]() {
			return y > bottomWall;
		}();
		particle.posY = hitsUp ? topWall : hitsDown ? bottomWall : y;
		particle.velY = hitsUp ? 0.0f : hitsDown ? 0.0f : particle.velY;
	};
	std::for_each(std::execution::par_unseq, particleMap.begin(), particleMap.end(), particleCollision);
}