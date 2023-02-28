module;

#include <algorithm>

export module Constants;

struct Border;

export namespace constants {
constexpr int maxwidth { 500 };
constexpr int maxheight { 500 };
constexpr double simheight { 500.0 };
constexpr double scale { (double)maxheight / simheight };
constexpr double simwidth { (double)maxwidth / scale };

enum class CellType
{
	Fluid,
	Air,
	Solid
};


//setupScene
constexpr double obstacleRadius = 20.00 * constants::scale;
constexpr double overRelaxation = 1.9;

constexpr double dt               = 1.0 / 60.0;
constexpr double numPressureIters = 50;
constexpr double numParticleIters = 2;

constexpr double gravity = -9.81 * scale;

constexpr double integ = dt * gravity;

constexpr double flipRatio           = 0.9;
constexpr double colorDiffusionCoeff = 0.001;
constexpr bool paused                = false;
constexpr bool showObstacle          = true;
constexpr bool showParticles         = true;
constexpr bool showGrid              = true;

constexpr double res = 50;

constexpr double tankHeight = simheight * scale;
constexpr double tankWidth  = simwidth * scale;
constexpr double h          = tankHeight / res;
constexpr double density    = 1000.0;

constexpr double relWaterHeight = 0.8;
constexpr double relWaterWidth  = 0.6;

// compute number of particles

constexpr double r  = 0.3 * h; // particle radius w.r.t. cell size
constexpr double dx = 2.0 * r;
const double dy     = sqrt(3.0) / 2.0 * dx;

const double numX         = floor((relWaterWidth * tankWidth - 2.0 * h - 2.0 * r) / dx);
const double numY         = floor((relWaterHeight * tankHeight - 2.0 * h - 2.0 * r) / dy);
const double maxParticles = numX * numY;

const int fNumX     = floor(tankWidth / h) + 1;
const int fNumY     = floor(tankHeight / h) + 1;
const int fNumCells = fNumX * fNumY;
// fluid
//
// particles
constexpr double pInvSpacing = 1.0 / (2.2 * r);
const int pNumX              = floor(tankWidth * pInvSpacing) + 1;
const int pNumY              = floor(tankHeight * pInvSpacing) + 1;
const int pNumCells          = pNumX * pNumY;

//Border fBorder { .maxX = fNumX, .maxY = fNumY };

} // namespace constants
  //double density, fInvSpacing, particleRestDensity, pInvSpacing, particleRadius, h, obstacleVelX, obstacleVelY, obstacleX, obstacleY;
  //int fNumX, fNumY, fNumCells, maxParticles, pNumX, pNumY, pNumCells;

//FlipFluid::FlipFluid(double density, double width, double height, double spacing, double particleRadius, int maxParticles) :
//   m_f = FlipFluid(density, tankWidth, tankHeight, h, r, maxParticles);
