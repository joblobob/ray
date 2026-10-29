module;

#include <algorithm>

export module Constants;

export namespace constants {
constexpr int maxwidth { 600 };
constexpr int maxheight { 400 };
constexpr float worldHeightMeters { 4.0f };

enum CellType { Fluid, Air, Solid, Soil, Rock };

constexpr int rockParticleCount = 1500;
constexpr int coarseSandParticleCount = 1500;
constexpr int fineSandParticleCount = 1500;
constexpr int grassParticleCount = 1500;

enum class ParticleMaterial { Rock, CoarseSand, FineSand, Grass, Water };

constexpr ParticleMaterial particleMaterialForId(const int id)
{
	if (id < rockParticleCount)
		return ParticleMaterial::Rock;
	if (id < rockParticleCount + coarseSandParticleCount)
		return ParticleMaterial::CoarseSand;
	if (id < rockParticleCount + coarseSandParticleCount + fineSandParticleCount)
		return ParticleMaterial::FineSand;
	if (id < rockParticleCount + coarseSandParticleCount + fineSandParticleCount + grassParticleCount)
		return ParticleMaterial::Grass;
	return ParticleMaterial::Water;
}

// Relative to water (1.0). Particle volume is unchanged, so these ratios
// also serve as relative particle masses in the grid transfer.
constexpr float particleDensityRatio(const ParticleMaterial material)
{
	switch (material) {
	case ParticleMaterial::Rock: return 3.50f;
	case ParticleMaterial::CoarseSand: return 2.66f;
	case ParticleMaterial::FineSand: return 2.00f;
	case ParticleMaterial::Grass: return 0.90f;
	case ParticleMaterial::Water: return 1.00f;
	}
	return 1.00f;
}

constexpr float particleDensityRatioForId(const int id)
{
	return particleDensityRatio(particleMaterialForId(id));
}

//setupScene

constexpr float dt = 1.0f / 60.0f;
constexpr int numPressureIters = 50;


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
constexpr float relWaterWidth  = 0.6f;

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
