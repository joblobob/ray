module;

#include <algorithm>
#include <string>
#include <unordered_map>
#include <vector>

#include <QString>

export module FlipFluid;

import Constants;
import CellCalculations;

export struct ExecutionLog {
	QString message;
	long long nsElapsed;
};

export struct Particle {
	double particlePosX = 0.0, particlePosY = 0.0, particleVelX = 0.0, particleVelY = 0.0;
	double particleColorR = 0.0, particleColorG = 0.0, particleColorB = 1.0;
};

export struct FlipFluid {
	double density, fInvSpacing, particleRestDensity, pInvSpacing, particleRadius, h, obstacleVelX, obstacleVelY, obstacleX, obstacleY;
	int fNumX, fNumY, fNumCells, maxParticles, pNumX, pNumY, pNumCells;
	std::vector<double> u, v, du, dv, prevU, prevV, s, cellColor, particleDensity;
	std::vector<int> numCellParticles, firstCellParticle, cellParticleIds;
	std::vector<constants::CellType> cellType;

	std::unordered_map<int, Particle> particleMap; // cell x particle

	Border pBorder;
	Border fBorder;

	FlipFluid() = default;
	FlipFluid(double density, double width, double height, double spacing, double particleRadius, int maxParticles);

	void setupObstacle(double x, double y, bool reset);
	void integrateParticles();
	void pushParticlesApart(int numIters);
	void handleParticleCollisions(double obstacleX, double obstacleY, double obstacleRadius);
	void updateParticleDensity();
	void transferVelocitiesToGrid();
	void transferVelocities(bool toGrid, double flipRatio);
	void solveIncompressibility(int numIters, double dt, double overRelaxation, bool compensateDrift = true);

	void updateParticleColors();
	void setSciColor(int cellNr, double val, double minVal, double maxVal);
	void updateCellColors();
	constexpr void calcColor(double& p1Color, double& p2Color);


	std::vector<ExecutionLog> simulate(double dt,
	    double flipRatio,
	    int numPressureIters,
	    int numParticleIters,
	    double overRelaxation,
	    bool compensateDrift,
	    bool separateParticles,
	    double obstacleRadius,
	    bool instrument);
};
