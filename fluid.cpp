#include "fluid.h"

#include <QtConcurrent>
#include <algorithm>
#include <execution>

#include <hitPosition.h>

#include <calc.h>

#include "ui_fluid.h"


void fluid::setupScene()
{
	float res = 100;

	float tankHeight = 1.0 * constants::simheight;
	float tankWidth  = 1.0 * constants::simwidth;
	float h          = tankHeight / res;
	float density    = 1000.0;

	float relWaterHeight = 0.8;
	float relWaterWidth  = 0.6;

	    // dam break

	    // compute number of particles

	float r = 0.3 * h; // particle radius w.r.t. cell size
	float dx    = 2.0 * r;
	float dy    = sqrt(3.0) / 2.0 * dx;

	float numX         = floor((relWaterWidth * tankWidth - 2.0 * h - 2.0 * r) / dx);
	float numY         = floor((relWaterHeight * tankHeight - 2.0 * h - 2.0 * r) / dy);
	float maxParticles = numX * numY;

	// create fluid

	f = FlipFluid(density, tankWidth, tankHeight, h, r, maxParticles);

	// create particles

	f.numParticles = numX * numY;
	int p = 0;
	for (int i = 0; i < numX; i++) {
		for (int j = 0; j < numY; j++) {
			f.particlePos[p++] = h + r + dx * i + (j % 2 == 0 ? 0.0 : r);
			f.particlePos[p++] = h + r + dy * j;
		}
	}

	// setup grid cells for tank

	float n = f.fNumY;

	for (int i = 0; i < f.fNumX; i++) {
		for (int j = 0; j < f.fNumY; j++) {
			float s = 1.0; // fluid
			if (i == 0 || i == f.fNumX - 1 || j == 0)
				s = 0.0; // solid
			f.s[i * n + j] = s;
		}
	}

	setupObstacle(3.0, 2.0, true);
}

void fluid::setupObstacle(int x, int y, bool reset)
{
	float vx = 0.0;
	float vy = 0.0;

	if (!reset) {
		vx = (x - obstacleX) / scenedt;
		vy = (y - obstacleY) / scenedt;
	}

	obstacleX = x;
	obstacleY = y;
	auto r           = constants::obstacleRadius;
	auto n           = f.fNumY;
	auto cd          = sqrt(2) * f.h;

	for (auto i = 1; i < f.fNumX - 2; i++) {
		for (auto j = 1; j < f.fNumY - 2; j++) {
			f.s[i * n + j] = 1.0;

			auto dx = (i + 0.5) * f.h - x;
			auto dy = (j + 0.5) * f.h - y;

			if (dx * dx + dy * dy < r * r) {
				f.s[i * n + j]       = 0.0;
				f.u[i * n + j]       = vx;
				f.u[(i + 1) * n + j] = vx;
				f.v[i * n + j]       = vy;
				f.v[i * n + j + 1]   = vy;
			}
		}
	}

	//scene.showObstacle = true;
	ggobstacleVelX = vx;
	ggobstacleVelY = vy;
}