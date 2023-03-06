module;

#include <algorithm>

export module Constants;

export namespace constants {
constexpr int maxwidth { 1200 };
constexpr int maxheight { 800 };

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
constexpr int numX         = static_cast<int>((relWaterWidth * maxwidth - 2.0 * cellHeight - 2.0 * particleRadius) / dx);
constexpr int numY         = static_cast<int>((relWaterHeight * maxheight - 2.0 * cellHeight - 2.0 * particleRadius) / dy);
constexpr int maxParticles = numX * numY;

constexpr int fNumX     = static_cast<int>(maxwidth / cellHeight);
constexpr int fNumY     = static_cast<int>(maxheight / cellHeight);
constexpr int fNumCells = fNumX * fNumY;
constexpr double fInvSpacing     = 1.0 / cellHeight;



} // namespace constants
