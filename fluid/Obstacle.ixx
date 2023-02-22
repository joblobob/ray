module;

#include <algorithm>
#include <execution>

export module Obstacle;

import BaseStructures;
import Constants;

export void setupObstacle(std::vector<Cell>& gridCells,
    const double newX,
    const double newY,
    double& obstacleX,
    double& obstacleY,
    double& obstacleVelX,
    double& obstacleVelY,
    const double h,
    const int fNumX,
    const int fNumY,
    bool reset)
{
	double vx = 0.0;
	double vy = 0.0;

	if (!reset) {
		vx = (newX - obstacleX) / constants::dt;
		vy = (newY - obstacleY) / constants::dt;
	}

	obstacleX = newX;
	obstacleY = newY;
	auto r    = constants::obstacleRadius;
	auto n    = fNumY;
	auto cd   = sqrt(2.0) * h;

	for (auto i = 1; i < fNumX - 2; i++) {
		for (auto j = 1; j < fNumY - 2; j++) {
			gridCells[i * n + j].s = 1.0;
			auto dx                = (i + 0.5) * h - newX;
			auto dy                = (j + 0.5) * h - newY;
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
	obstacleVelX = vx;
	obstacleVelY = vy;
}