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

constexpr float res = 100.0f;

constexpr float cellHeight     = maxheight / res;
constexpr float halfCellHeight = 0.5f * constants::cellHeight;

constexpr float density = 1000.0f;



// Particle layout is in scene pixels; changing grid resolution does not change particle size.
constexpr float particleRadius       = 1.2f;
constexpr float particleDiameter     = 2.0f * particleRadius;
constexpr float particleSpacingX    = particleDiameter;
constexpr float particleRowSpacing  = 1.7320508f * particleRadius;
constexpr float particleStaggerX    = particleRadius;

constexpr float relWaterHeight = 0.8f;
constexpr float relWaterWidth  = 0.8f;

constexpr float particleStartX = cellHeight + particleRadius;
constexpr float particleStartY = cellHeight + particleRadius;
constexpr float availableParticleWidth = relWaterWidth * maxwidth - 2.0f * particleStartX;
constexpr float availableParticleHeight = relWaterHeight * maxheight - 2.0f * particleStartY;
constexpr int numX = static_cast<int>(availableParticleWidth / particleSpacingX);
constexpr int numY = static_cast<int>(availableParticleHeight / particleRowSpacing);
constexpr int maxParticles     = numX * numY;

constexpr int fNumX         = static_cast<int>(maxwidth / cellHeight);
constexpr int fNumY         = static_cast<int>(maxheight / cellHeight);
constexpr int fNumCells     = fNumX * fNumY;
constexpr float fInvSpacing = 1.0f / cellHeight;



} // namespace constants
