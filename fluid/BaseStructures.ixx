module;

#include <algorithm>
#include <cstddef>
#include <vector>

export module BaseStructures;

import Constants;

export struct Border {
	int minX = 0, minY = 0, maxX = constants::maxwidth, maxY = constants::maxwidth;
};

// Hot particle fields live in separate contiguous arrays. Array indices are
// stable particle IDs used by the spatial grid and particle buffers.
export struct ParticleStore {
	std::vector<float> posX;
	std::vector<float> posY;
	std::vector<float> velX;
	std::vector<float> velY;
	std::vector<float> colorR;
	std::vector<float> colorG;
	std::vector<float> colorB;

	explicit ParticleStore(const std::size_t count = 0) { resize(count); }

	void resize(const std::size_t count)
	{
		posX.resize(count);
		posY.resize(count);
		velX.resize(count);
		velY.resize(count);
		colorR.resize(count);
		colorG.resize(count);
		colorB.resize(count);
	}

	[[nodiscard]] std::size_t size() const { return posX.size(); }
};


export struct Cell {
	int cellNumX, cellNumY;
	float u, v, du, dv, prevU, prevV, s, particleDensity = 0.0f, particleMassDensity = 0.0f;
	float colorR = 0.0f, colorG = 0.0f, colorB = 0.0f;
	constants::CellType cellType = constants::CellType::Solid;
};
