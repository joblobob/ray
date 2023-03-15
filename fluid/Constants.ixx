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

constexpr float dt = 1.0f / 60.0f;


constexpr bool paused        = false;
constexpr bool showObstacle  = true;
constexpr bool showParticles = true;
constexpr bool showGrid      = false;

constexpr float res = 120.0;

constexpr float cellHeight     = maxheight / res;
constexpr float halfCellHeight = 0.5f * constants::cellHeight;

constexpr float density = 1000.0f;



// compute number of particles

constexpr float particleRadius = 0.3f * cellHeight; // particle radius w.r.t. cell size
constexpr float dx             = 2.0f * particleRadius;
constexpr float sqrt_of_3 { 1.73205080f };
constexpr float dy = sqrt_of_3 / 2.0f * dx;

constexpr float relWaterHeight = 0.8f;
constexpr float relWaterWidth  = 0.6f;
constexpr int numX             = static_cast<int>((relWaterWidth * maxwidth - 2.0f * cellHeight - 2.0f * particleRadius) / dx);
constexpr int numY             = static_cast<int>((relWaterHeight * maxheight - 2.0f * cellHeight - 2.0f * particleRadius) / dy);
constexpr int maxParticles     = numX * numY;

constexpr int fNumX         = static_cast<int>(maxwidth / cellHeight);
constexpr int fNumY         = static_cast<int>(maxheight / cellHeight);
constexpr int fNumCells     = fNumX * fNumY;
constexpr float fInvSpacing = 1.0f / cellHeight;



} // namespace constants
