module;

#include <algorithm>
#include <vector>

export module FlipFluid;

export namespace constants {
constexpr int maxwidth { 500 };
constexpr int maxheight { 500 };
constexpr double simheight { 3.0 };
inline double scale { (double)maxheight / simheight };
inline double simwidth { (double)maxwidth / scale };

constexpr int U_FIELD = 0;
constexpr int V_FIELD = 1;

constexpr int FLUID_CELL = 0;
constexpr int AIR_CELL   = 1;
constexpr int SOLID_CELL = 2;

constexpr int cnt = 0;

//setupScene
const double obstacleRadius     = 0.20 * constants::scale;
constexpr double overRelaxation = 1.9;

constexpr double dt               = 1.0 / 60.0;
constexpr double numPressureIters = 100;
constexpr double numParticleIters = 2;

const double gravity = -9.81 * scale;
//constexpr double dt : 1.0 / 120.0,
constexpr double flipRatio       = 0.9;
constexpr bool compensateDrift   = true;
constexpr bool separateParticles = true;
constexpr bool paused            = false;
constexpr bool showObstacle      = true;
constexpr bool showParticles     = true;
constexpr bool showGrid          = true;

} // namespace constants



export struct FlipFluid {
	double density, fInvSpacing, particleRestDensity, pInvSpacing, particleRadius, h, obstacleVelX, obstacleVelY, obstacleX, obstacleY;
	int fNumX, fNumY, fNumCells, maxParticles, pNumX, pNumY, pNumCells, numParticles;
	std::vector<double> u, v, du, dv, prevU, prevV, s, cellColor, particlePos, particleColor, particleVel, particleDensity;
	std::vector<int> cellType, numCellParticles, firstCellParticle, cellParticleIds;

	FlipFluid() = default;
	FlipFluid(double density, double width, double height, double spacing, double particleRadius, int maxParticles);
	inline bool isVeryCloseToZero(double x);

	constexpr int cellNumberFromParticle(double x, double y, int clampMin = 0);
	constexpr int cellNumberFromGrid(double x, double y, int clampMin = 0);

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
};
