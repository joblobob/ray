export module Constants;

export namespace constants {
constexpr int maxwidth { 500 };
constexpr int maxheight { 500 };
constexpr double simheight { 3.0 };
constexpr double scale { (double)maxheight / simheight };
constexpr double simwidth { (double)maxwidth / scale };

enum class CellType
{
	Fluid,
	Air,
	Solid
};

//setupScene
constexpr double obstacleRadius = 0.20 * constants::scale;
constexpr double overRelaxation = 1.9;

constexpr double dt               = 1.0 / 60.0;
constexpr double numPressureIters = 50;
constexpr double numParticleIters = 2;

const double gravity = -9.81 * scale;

constexpr double flipRatio       = 0.9;
constexpr double colorDiffusionCoeff    = 0.001;
constexpr bool compensateDrift   = true;
constexpr bool separateParticles = true;
constexpr bool paused            = false;
constexpr bool showObstacle      = true;
constexpr bool showParticles     = true;
constexpr bool showGrid          = true;

} // namespace constants