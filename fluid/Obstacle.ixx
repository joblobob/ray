module;

#include <algorithm>
#include <execution>

export module Obstacle;

import BaseStructures;
import Constants;

export namespace ObstacleConstants {
constexpr double radius { 20.00 };
}

export struct Obstacle {
	double x = 0.0, y = 0.0, velX = 0.0, velY = 0.0;
};

export void setupObstacle(std::vector<Cell>& gridCells, const double newX, const double newY, Obstacle& obstacle, bool reset)
{
	double vx = 0.0;
	double vy = 0.0;

	if (!reset) {
		vx = (newX - obstacle.x) / constants::dt;
		vy = (newY - obstacle.y) / constants::dt;
	}

	obstacle.x           = newX;
	obstacle.y           = newY;
	const double r       = ObstacleConstants::radius;
	const int n = constants::fNumY;

	for (int i = 1; i < constants::fNumX - 2; i++) {
		for (int j = 1; j < constants::fNumY - 2; j++) {
			gridCells[i * n + j].s = 1.0;
			const double dx        = (i + 0.5) * constants::cellHeight - newX;
			const double dy        = (j + 0.5) * constants::cellHeight - newY;
			if (dx * dx + dy * dy < r * r) {
				gridCells[i * n + j].s       = 0.0;
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