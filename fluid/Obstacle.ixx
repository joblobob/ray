module;

#include <algorithm>
#include <execution>

export module Obstacle;

import BaseStructures;
import Constants;

export namespace ObstacleConstants {
constexpr float radius { 20.00f };
}

export struct Obstacle {
	float x = 0.0f, y = 0.0f, velX = 0.0f, velY = 0.0f;
};

export void setupObstacle(std::vector<Cell>& gridCells, const float newX, const float newY, Obstacle& obstacle, bool reset)
{
	float vx = 0.0f;
	float vy = 0.0f;

	if (!reset) {
		vx = (newX - obstacle.x) / constants::dt;
		vy = (newY - obstacle.y) / constants::dt;
	}

	obstacle.x    = newX;
	obstacle.y    = newY;
	const float r = ObstacleConstants::radius;
	const int n   = constants::fNumY;

	for (int i = 1; i < constants::fNumX - 2; i++) {
		for (int j = 1; j < constants::fNumY - 2; j++) {
			gridCells[i * n + j].s = 1.0f;
			const float dx         = (i + 0.5f) * constants::cellHeight - newX;
			const float dy         = (j + 0.5f) * constants::cellHeight - newY;
			if (dx * dx + dy * dy < r * r) {
				gridCells[i * n + j].s       = 0.0f;
				gridCells[i * n + j].u       = vx;
				gridCells[(i + 1) * n + j].u = vx;
				gridCells[i * n + j].v       = vy;
				gridCells[i * n + j + 1].v   = vy;
			}
		}
	}

	//scene.showObstacle = true;
	obstacle.velX = vx;
	obstacle.velY = vy;
}