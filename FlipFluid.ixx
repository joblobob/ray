module;

#include <algorithm>
#include <string>
#include <unordered_map>
#include <vector>

#include <QString>

export module FlipFluid;

import BaseStructures;
import Constants;
import CellCalculations;

export struct ExecutionLog {
	QString message;
	long long nsElapsed;
};

export struct FlipFluid {
	double density, fInvSpacing, particleRestDensity, pInvSpacing, particleRadius, h, obstacleVelX, obstacleVelY, obstacleX, obstacleY;
	int fNumX, fNumY, fNumCells, maxParticles, pNumX, pNumY, pNumCells;
	std::vector<int> cellParticleIds;

	std::vector<Particle> particleMap; // particles

	std::vector<ParticleInCells> particleCells; // particlesincellinfo

	std::vector<Cell> gridCells;

	Border pBorder;
	Border fBorder;

	FlipFluid() = default;
	FlipFluid(double density, double width, double height, double spacing, double particleRadius, int maxParticles);

	void setupObstacle(double x, double y, bool reset);
	void pushParticlesApart(int numIters);
	void handleParticleCollisions(double obstacleX, double obstacleY, double obstacleRadius);
	void updateParticleDensity();
	void transferVelocitiesToGrid();
	void transferVelocitiesToParticles();
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

private:
	void parseVelocitiesV(Particle& particle);
	void restoreSolidCellsV(Cell& cell);

	void parseVelocitiesParticles(Particle& particle);

	void setVelComponent(double& f, double& d, double pv, double delta);
};
