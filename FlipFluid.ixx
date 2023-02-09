module;

#include <algorithm>
#include <string>
#include <vector>

#include <QString>

export module FlipFluid;

import Constants;
import CellCalculations;

export struct ExecutionLog {
	QString message;
	long long nsElapsed;
};


export struct FlipFluid {
	double density, fInvSpacing, particleRestDensity, pInvSpacing, particleRadius, h, obstacleVelX, obstacleVelY, obstacleX, obstacleY;
	int fNumX, fNumY, fNumCells, maxParticles, pNumX, pNumY, pNumCells, numParticles;
	std::vector<double> u, v, du, dv, prevU, prevV, s, cellColor, particlePos, particleColor, particleVel, particleDensity;
	std::vector<int> numCellParticles, firstCellParticle, cellParticleIds;
	std::vector<constants::CellType> cellType;

	FlipFluid() = default;
	FlipFluid(double density, double width, double height, double spacing, double particleRadius, int maxParticles);

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


	std::vector<ExecutionLog> simulate(double dt,
	    double gravity,
	    double flipRatio,
	    int numPressureIters,
	    int numParticleIters,
	    double overRelaxation,
	    bool compensateDrift,
	    bool separateParticles,
	    double obstacleRadius,
	    bool instrument);

private:
	Border pBorder;
	Border fBorder;
};
