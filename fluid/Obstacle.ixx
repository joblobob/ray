module;

#include "../mdspan.hpp"
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

	std::experimental::mdspan gridView(gridCells.data(), constants::fNumX, constants::fNumY);

	for (std::size_t i = 1; i < gridView.extent(0) - 2; i++) {
		for (std::size_t j = 1; j < gridView.extent(1) - 2; j++) {
			gridView(i, j).s = 1.0f;
			const float dx   = (i + 0.5f) * constants::cellHeight - newX;
			const float dy   = (j + 0.5f) * constants::cellHeight - newY;
			if (dx * dx + dy * dy < r * r) {
				gridView(i, j).s       = 0.0f;
				gridView(i, j).u       = vx;
				gridView((i + 1), j).u = vx;
				gridView(i, j).v       = vy;
				gridView(i, j + 1).v   = vy;
			}
		}
	}


	obstacle.velX = vx;
	obstacle.velY = vy;
}