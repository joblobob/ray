module;

#include <algorithm>
#include <vector>

export module FlipFluid;

import Constants;
import CellCalculations;


export struct FlipFluid {
	double density, fInvSpacing, particleRestDensity, pInvSpacing, particleRadius, h, obstacleVelX, obstacleVelY, obstacleX, obstacleY;
	int fNumX, fNumY, fNumCells, maxParticles, pNumX, pNumY, pNumCells, numParticles;
	std::vector<double> u, v, du, dv, prevU, prevV, s, cellColor, particlePos, particleColor, particleVel, particleDensity;
	std::vector<int> cellType, numCellParticles, firstCellParticle, cellParticleIds;

	FlipFluid() = default;
	FlipFluid(double density, double width, double height, double spacing, double particleRadius, int maxParticles);
	inline bool isVeryCloseToZero(double x);

	void setupObstacle(double x, double y, bool reset);
	void integrateParticles(double dt, double gravity);
	void pushParticlesApart(int numIters);
	void handleParticleCollisions(double obstacleX, double obstacleY, double obstacleRadius);
	void updateParticleDensity();
	void transferVelocities(bool toGrid, double flipRatio);
	void solveIncompressibility(int numIters, double dt, double overRelaxation, bool compensateDrift = true);

	void updateParticleColors();
	void setSciColor(int cellNr, double val, double minVal, double maxVal);
	void updateCellColors();
	

	void simulate(double dt, double gravity, double flipRatio, int numPressureIters, int numParticleIters, double overRelaxation,
	    bool compensateDrift,
	    bool separateParticles,
	    double obstacleRadius);

	private:
		Border pBorder;
		Border fBorder;
};
