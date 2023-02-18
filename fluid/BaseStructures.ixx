export module BaseStructures;

import Constants;

export struct Particle {
	int id;
	double posX = 0.0, posY = 0.0, velX = 0.0, velY = 0.0;
	double colorR = 0.0, colorG = 0.0, colorB = 1.0;
};

export struct ParticleInCells {
	int numCellParticles = 0, firstCellParticle = 0;
};

export struct Cell {
	int cellNumX, cellNumY;
	double u, v, du, dv, prevU, prevV, s, particleDensity = 0.0;
	double colorR = 0.0, colorG = 0.0, colorB = 0.0;
	constants::CellType cellType = constants::CellType::Solid;
};