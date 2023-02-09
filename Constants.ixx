export module Constants;

export namespace constants
{
	constexpr int maxwidth { 500 };
	constexpr int maxheight { 500 };
	constexpr double simheight { 3.0 };
	inline double scale { (double)maxheight / simheight };
	inline double simwidth { (double)maxwidth / scale };

	constexpr int U_FIELD = 0;
	constexpr int V_FIELD = 1;

	constexpr int FLUID_CELL = 0;
	constexpr int AIR_CELL   = 1;
	constexpr int SOLID_CELL = 2;

	constexpr int cnt = 0;

	//setupScene
	const double obstacleRadius     = 0.20 * constants::scale;
	constexpr double overRelaxation = 1.9;

	constexpr double dt               = 1.0 / 60.0;
	constexpr double numPressureIters = 100;
	constexpr double numParticleIters = 2;

	const double gravity = -9.81 * scale;
	//constexpr double dt : 1.0 / 120.0,
	constexpr double flipRatio       = 0.9;
	constexpr bool compensateDrift   = true;
	constexpr bool separateParticles = true;
	constexpr bool paused            = false;
	constexpr bool showObstacle      = true;
	constexpr bool showParticles     = true;
	constexpr bool showGrid          = true;

} // namespace constants