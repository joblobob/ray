export module BaseStructures;

import Constants;

export struct Border {
	int minX = 0, minY = 0, maxX = constants::maxwidth, maxY = constants::maxwidth;
};

export struct Particle {
	int id;
	float posX = 0.0f, posY = 0.0f, velX = 0.0f, velY = 0.0f;
	float colorR = 0.0f, colorG = 0.0f, colorB = 1.0f;
};


export struct Cell {
	int cellNumX, cellNumY;
	float u, v, du, dv, prevU, prevV, s, particleDensity = 0.0f;
	float colorR = 0.0f, colorG = 0.0f, colorB = 0.0f;
	constants::CellType cellType = constants::CellType::Solid;
};