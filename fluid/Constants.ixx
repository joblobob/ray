module;

#include <algorithm>

export module Constants;

export namespace constants {
constexpr unsigned int maxwidth { 800 };
constexpr unsigned int maxheight { 500 };

enum class CellType
{
	Fluid,
	Air,
	Solid
};

//setupScene

constexpr double dt = 1.0 / 20.0;


constexpr bool paused        = false;
constexpr bool showObstacle  = true;
constexpr bool showParticles = true;
constexpr bool showGrid      = false;

constexpr double res = 100;

constexpr double cellHeight     = maxheight / res;
constexpr double halfCellHeight = 0.5 * constants::cellHeight;

constexpr double density = 1000.0;



// compute number of particles

constexpr double particleRadius = 0.3 * cellHeight; // particle radius w.r.t. cell size
constexpr double dx             = 2.0 * particleRadius;
constexpr double sqrt_of_3 { 1.7320508075688772 };
constexpr double dy = sqrt_of_3 / 2.0 * dx;

constexpr double relWaterHeight     = 0.8;
constexpr double relWaterWidth      = 0.6;
constexpr unsigned int numX         = static_cast<unsigned int>((relWaterWidth * maxwidth - 2.0 * cellHeight - 2.0 * particleRadius) / dx);
constexpr unsigned int numY         = static_cast<unsigned int>((relWaterHeight * maxheight - 2.0 * cellHeight - 2.0 * particleRadius) / dy);
constexpr unsigned int maxParticles = numX * numY;

constexpr unsigned int fNumX     = static_cast<int>(maxwidth / cellHeight);
constexpr unsigned int fNumY     = static_cast<int>(maxheight / cellHeight);
constexpr unsigned int fNumCells = fNumX * fNumY;
constexpr double fInvSpacing     = 1.0 / cellHeight;



} // namespace constants
