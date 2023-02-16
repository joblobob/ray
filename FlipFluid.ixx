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
	int id;
	double posX = 0.0, posY = 0.0, velX = 0.0, velY = 0.0;
	double colorR = 0.0, colorG = 0.0, colorB = 1.0;
};

struct ParticleInCells {
	int numCellParticles = 0, firstCellParticle = 0;
};

export struct Cell {
	int cellNumX, cellNumY;
	double u, v, du, dv, prevU, prevV, s, particleDensity = 0.0;
	double colorR = 0.0, colorG = 0.0, colorB = 0.0;
	constants::CellType cellType = constants::CellType::Solid;
};

export struct FlipFluid {
	double density, fInvSpacing, particleRestDensity, pInvSpacing, particleRadius, h, obstacleVelX, obstacleVelY, obstacleX, obstacleY;
	int fNumX, fNumY, fNumCells, maxParticles, pNumX, pNumY, pNumCells;
	//std::vector<double> u, v, du, dv, prevU, prevV, s, cellColor, particleDensity;
	std::vector<int> cellParticleIds;
	//std::vector<constants::CellType> cellType;

	std::vector<Particle> particleMap; // particles

	std::vector<ParticleInCells> particleCells; // particlesincellinfo

	std::vector<Cell> gridCells;

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
	void setSciColor(Cell& cell, double val, double minVal, double maxVal);
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
