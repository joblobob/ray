module;

#include <algorithm>

export module Constants;

export namespace constants {
constexpr int maxwidth { 600 };
constexpr int maxheight { 400 };
constexpr float worldHeightMeters { 4.0f };

enum CellType { Fluid, Air, Solid, Soil, Rock };

//setupScene

constexpr float dt = 1.0f / 60.0f;


constexpr bool paused        = false;
constexpr bool showObstacle  = true;
constexpr bool showParticles = true;
constexpr bool showGrid      = false;

constexpr float res = 50.0f;

constexpr float cellHeight     = maxheight / res;
constexpr float halfCellHeight = 0.5f * constants::cellHeight;

constexpr float density = 1000.0f;



// compute number of particles

// Keep the visible particle size and spacing independent of the simulation-grid resolution.
constexpr float particleSizeReferenceRes = 100.0f;
constexpr float particleRadius = 0.3f * (maxheight / particleSizeReferenceRes);
constexpr float dx             = 2.0f * particleRadius;
constexpr float sqrt_of_3 { 1.73205080f };
constexpr float dy = sqrt_of_3 / 2.0f * dx;

// Seed a larger volume without changing particle spacing or radius.
constexpr float relWaterHeight = 1.f;
constexpr float relWaterWidth  = 1.0f;
constexpr int numX             = static_cast<int>((relWaterWidth * maxwidth - 2.0f * cellHeight - 2.0f * particleRadius) / dx);
constexpr int numY             = static_cast<int>((relWaterHeight * maxheight - 2.0f * cellHeight - 2.0f * particleRadius) / dy);
constexpr int maxParticles     = numX * numY;

constexpr int fNumX         = static_cast<int>(maxwidth / cellHeight);
constexpr int fNumY         = static_cast<int>(maxheight / cellHeight);
constexpr int fNumCells     = fNumX * fNumY;
constexpr float fInvSpacing = 1.0f / cellHeight;



} // namespace constants
