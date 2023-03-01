module;

#include <algorithm>

export module Constants;

struct Border;

export constexpr int int_floor(const double d)
{
	const int i = static_cast<int>(d);
	return d < i ? i - 1 : i;
}

export namespace constants {
constexpr int maxwidth { 500 };
constexpr int maxheight { 500 };

enum class CellType
{
	Fluid,
	Air,
	Solid
};


//setupScene
constexpr double obstacleRadius = 20.00;
constexpr double overRelaxation = 1.9;

constexpr double dt               = 1.0 / 5.0;
constexpr double numPressureIters = 50;
constexpr double numParticleIters = 2;

constexpr double gravity = -9.81;

constexpr double integ = dt * gravity;

constexpr double flipRatio           = 0.9;
constexpr double colorDiffusionCoeff = 0.001;
constexpr bool paused                = false;
constexpr bool showObstacle          = true;
constexpr bool showParticles         = true;
constexpr bool showGrid              = false;

constexpr double res = 100;

constexpr double h       = maxheight / res;
constexpr double density = 1000.0;

constexpr double relWaterHeight = 0.8;
constexpr double relWaterWidth  = 0.6;

// compute number of particles

constexpr double r  = 0.3 * h; // particle radius w.r.t. cell size
constexpr double dx = 2.0 * r;
constexpr double sqrt_of_3 { 1.7320508075688772 };
constexpr double sqrt_of_2 { 1.4142135623730950 };
constexpr double dy = sqrt_of_3 / 2.0 * dx;

constexpr double numX         = int_floor((relWaterWidth * maxwidth - 2.0 * h - 2.0 * r) / dx);
constexpr double numY         = int_floor((relWaterHeight * maxheight - 2.0 * h - 2.0 * r) / dy);
constexpr double maxParticles = numX * numY;

constexpr int fNumX     = int_floor(maxwidth / h) + 1;
constexpr int fNumY     = int_floor(maxheight / h) + 1;
constexpr int fNumCells = fNumX * fNumY;
// fluid
//
// particles
constexpr double pInvSpacing = 1.0 / (2.2 * r);
constexpr int pNumX          = int_floor(maxwidth * pInvSpacing) + 1;
constexpr int pNumY          = int_floor(maxheight * pInvSpacing) + 1;
constexpr int pNumCells      = pNumX * pNumY;


} // namespace constants
