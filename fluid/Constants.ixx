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

constexpr double dt = 1.0 / 5.0;


constexpr bool paused        = false;
constexpr bool showObstacle  = true;
constexpr bool showParticles = true;
constexpr bool showGrid      = false;

constexpr double res = 100;

constexpr double cellHeight     = maxheight / res;
constexpr double halfCellHeight = 0.5 * constants::cellHeight;

constexpr double density = 1000.0;
constexpr double k       = 1.0;

constexpr double relWaterHeight = 0.8;
constexpr double relWaterWidth  = 0.6;

// compute number of particles

constexpr double particleRadius = 0.3 * cellHeight; // particle radius w.r.t. cell size
constexpr double dx             = 2.0 * particleRadius;
constexpr double sqrt_of_3 { 1.7320508075688772 };
constexpr double sqrt_of_2 { 1.4142135623730950 };
constexpr double dy = sqrt_of_3 / 2.0 * dx;

constexpr double numX         = int_floor((relWaterWidth * maxwidth - 2.0 * cellHeight - 2.0 * particleRadius) / dx);
constexpr double numY         = int_floor((relWaterHeight * maxheight - 2.0 * cellHeight - 2.0 * particleRadius) / dy);
constexpr double maxParticles = numX * numY;

constexpr int fNumX          = int_floor(maxwidth / cellHeight);
constexpr int fNumY          = int_floor(maxheight / cellHeight);
constexpr int fNumCells      = fNumX * fNumY;
constexpr double fInvSpacing = 1.0 / cellHeight;
// fluid
//
// particles
constexpr double pInvSpacing = 1.0 / (2.2 * particleRadius);
constexpr int pNumX          = int_floor(maxwidth * pInvSpacing) + 1;
constexpr int pNumY          = int_floor(maxheight * pInvSpacing) + 1;
constexpr int pNumCells      = pNumX * pNumY;


} // namespace constants
