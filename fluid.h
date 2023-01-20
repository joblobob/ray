#ifndef FLUID_H
#define FLUID_H

#include "qelapsedtimer.h"

#include <QDialog>
#include <QFile>
#include <QGraphicsPixmapItem>
#include <QGraphicsScene>
#include <QImage>


namespace constants {
constexpr int maxwidth { 300 };
//constexpr int maxheight { 200 };
constexpr float simheight { 300 };
inline float scale { 1 };
inline float simwidth { (float)maxwidth / scale };

constexpr int U_FIELD = 0;
constexpr int V_FIELD = 1;

constexpr int FLUID_CELL = 0;
constexpr int AIR_CELL   = 1;
constexpr int SOLID_CELL = 2;

constexpr int cnt = 0;

//setupScene
constexpr float obstacleRadius = 50;
constexpr float overRelaxation = 1.9;

constexpr float dt               = 1.0 / 60.0f;
constexpr float numPressureIters = 50;
constexpr float numParticleIters = 2;

constexpr float gravity = -9.81;
//constexpr float dt : 1.0 / 120.0,
constexpr float flipRatio        = 0.9;
constexpr bool compensateDrift   = true;
constexpr bool separateParticles = true;
constexpr bool paused            = false;
constexpr bool showObstacle      = true;
constexpr float obstacleVelX     = 0.0;
constexpr float obstacleVelY     = 0.0;
constexpr bool showParticles     = true;
constexpr bool showGrid          = true;

} // namespace constants

namespace {
float ggobstacleVelX = 0.0;
float ggobstacleVelY = 0.0;
} // namespace


struct FlipFluid {
	float density, fInvSpacing, particleRestDensity, pInvSpacing, particleRadius, h;
	int fNumX, fNumY, fNumCells, maxParticles, pNumX, pNumY, pNumCells, numParticles;
	std::vector<float> u, v, du, dv, prevU, prevV, p, s, cellColor, particlePos, particleColor, particleVel, particleDensity;
	std::vector<int> cellType, numCellParticles, firstCellParticle, cellParticleIds;

	FlipFluid() = default;
	FlipFluid(float density, float width, float height, float spacing, float particleRadius, int maxParticles) :
	    density(density),
	    fNumX(floor(width / spacing) + 1),
	    fNumY(floor(height / spacing) + 1),
	    fNumCells(fNumX * fNumY),
	    u(fNumCells),
	    // fluid
	    v(fNumCells),
	    du(fNumCells),
	    dv(fNumCells),
	    prevU(fNumCells),
	    prevV(fNumCells),
	    p(fNumCells),
	    s(fNumCells),
	    cellType(fNumCells),
	    cellColor(3 * fNumCells),
	    particleRadius(particleRadius),
	    // particles
	    maxParticles(maxParticles),
	    particlePos(2 * maxParticles),
	    particleColor(3 * maxParticles),
	    particleVel(2 * maxParticles),
	    particleDensity(fNumCells),
	    particleRestDensity(0.0),
	    pInvSpacing(1.0 / (2.2 * particleRadius)),
	    pNumX(floor(width * pInvSpacing) + 1),
	    pNumY(floor(height * pInvSpacing) + 1),
	    pNumCells(pNumX * pNumY),
	    numCellParticles(pNumCells),
	    firstCellParticle(pNumCells + 1),
	    cellParticleIds(maxParticles),
	    numParticles(0)
	{
		h           = std::max(width / (float)fNumX, height / (float)fNumY);
		fInvSpacing = 1.0 / h;
		//set default color to blue
		for (auto i = 0; i < maxParticles; i++)
			particleColor[3 * i + 2] = 1.0;
	}

	void integrateParticles(float dt, float gravity)
	{
		for (auto i = 0; i < numParticles; i++) {
			particleVel[2 * i + 1] += dt * gravity;
			particlePos[2 * i] += particleVel[2 * i] * dt;
			particlePos[2 * i + 1] += particleVel[2 * i + 1] * dt;
		}
	}

	void pushParticlesApart(int numIters)
	{
		auto colorDiffusionCoeff = 0.001f;

		// count particles per cell

		std::fill(numCellParticles.begin(), numCellParticles.end(), 0);

		for (auto i = 0; i < numParticles; i++) {
			auto x = particlePos[2 * i];
			auto y = particlePos[2 * i + 1];

			auto xi     = std::clamp(floor(x * pInvSpacing), 0.0f, (float)pNumX - 1);
			auto yi     = std::clamp(floor(y * pInvSpacing), 0.0f, (float)pNumY - 1);
			auto cellNr = xi * pNumY + yi;
			numCellParticles[cellNr]++;
		}

		// partial sums

		auto first = 0;

		for (auto i = 0; i < pNumCells; i++) {
			first += numCellParticles[i];
			firstCellParticle[i] = first;
		}
		firstCellParticle[pNumCells] = first; // guard

		// fill particles into cells

		for (auto i = 0; i < numParticles; i++) {
			auto x = particlePos[2 * i];
			auto y = particlePos[2 * i + 1];

			auto xi     = std::clamp(floor(x * pInvSpacing), 0.0f, (float)pNumX - 1);
			auto yi     = std::clamp(floor(y * pInvSpacing), 0.0f, (float)pNumY - 1);
			auto cellNr = xi * pNumY + yi;
			firstCellParticle[cellNr]--;
			cellParticleIds[firstCellParticle[cellNr]] = i;
		}

		// push particles apart

		auto minDist  = 2.0 * particleRadius;
		auto minDist2 = minDist * minDist;

		for (auto iter = 0; iter < numIters; iter++) {
			for (auto i = 0; i < numParticles; i++) {
				auto px = particlePos[2 * i];
				auto py = particlePos[2 * i + 1];

				int pxi = floor(px * pInvSpacing);
				int pyi = floor(py * pInvSpacing);
				auto x0 = std::max(pxi - 1, 0);
				auto y0 = std::max(pyi - 1, 0);
				auto x1 = std::min(pxi + 1, pNumX - 1);
				auto y1 = std::min(pyi + 1, pNumY - 1);

				for (auto xi = x0; xi <= x1; xi++) {
					for (auto yi = y0; yi <= y1; yi++) {
						auto cellNr = xi * pNumY + yi;
						auto first  = firstCellParticle[cellNr];
						auto last   = firstCellParticle[cellNr + 1];
						for (auto j = first; j < last; j++) {
							auto id = cellParticleIds[j];
							if (id == i)
								continue;
							auto qx = particlePos[2 * id];
							auto qy = particlePos[2 * id + 1];

							auto dx = qx - px;
							auto dy = qy - py;
							auto d2 = dx * dx + dy * dy;
							if (d2 > minDist2 || d2 == 0.0)
								continue;
							auto d = sqrt(d2);
							auto s = 0.5 * (minDist - d) / d;
							dx *= s;
							dy *= s;
							particlePos[2 * i] -= dx;
							particlePos[2 * i + 1] -= dy;
							particlePos[2 * id] += dx;
							particlePos[2 * id + 1] += dy;

							// diffuse colors

							for (auto k = 0; k < 3; k++) {
								auto color0               = particleColor[3 * i + k];
								auto color1               = particleColor[3 * id + k];
								auto color                = (color0 + color1) * 0.5;
								particleColor[3 * i + k]  = color0 + (color - color0) * colorDiffusionCoeff;
								particleColor[3 * id + k] = color1 + (color - color1) * colorDiffusionCoeff;
							}
						}
					}
				}
			}
		}
	}

	void handleParticleCollisions(int obstacleX, int obstacleY, float obstacleRadius)
	{
		auto h        = 1.0 / fInvSpacing;
		auto r        = particleRadius;
		auto or2      = obstacleRadius * obstacleRadius;
		auto minDist  = obstacleRadius + r;
		auto minDist2 = minDist * minDist;

		auto minX = h + r;
		auto maxX = (fNumX - 1) * h - r;
		auto minY = h + r;
		auto maxY = (fNumY - 1) * h - r;


		for (auto i = 0; i < numParticles; i++) {
			auto x = particlePos[2 * i];
			auto y = particlePos[2 * i + 1];

			auto dx = x - obstacleX;
			auto dy = y - obstacleY;
			auto d2 = dx * dx + dy * dy;

			// obstacle collision

			if (d2 < minDist2) {
				// auto d = Math.sqrt(d2);
				// auto s = (minDist - d) / d;
				// x += dx * s;
				// y += dy * s;

				particleVel[2 * i]     = ggobstacleVelX;
				particleVel[2 * i + 1] = ggobstacleVelY;
			}

			// wall collisions

			if (x < minX) {
				x                  = minX;
				particleVel[2 * i] = 0.0;
			}
			if (x > maxX) {
				x                  = maxX;
				particleVel[2 * i] = 0.0;
			}
			if (y < minY) {
				y                      = minY;
				particleVel[2 * i + 1] = 0.0;
			}
			if (y > maxY) {
				y                      = maxY;
				particleVel[2 * i + 1] = 0.0;
			}
			particlePos[2 * i]     = x;
			particlePos[2 * i + 1] = y;
		}
	}

	void updateParticleDensity()
	{
		auto n   = fNumY;
		auto h1  = fInvSpacing;
		float h2 = 0.5 * h;


		std::fill(particleDensity.begin(), particleDensity.end(), 0.0f);


		for (auto i = 0; i < numParticles; i++) {
			auto x = particlePos[2 * i];
			auto y = particlePos[2 * i + 1];

			x = std::clamp(x, h, (fNumX - 1) * h);
			y = std::clamp(y, h, (fNumY - 1) * h);

			auto x0 = floor((x - h2) * h1);
			auto tx = ((x - h2) - x0 * h) * h1;
			auto x1 = std::min(x0 + 1, (float)fNumX - 2);

			auto y0 = floor((y - h2) * h1);
			auto ty = ((y - h2) - y0 * h) * h1;
			auto y1 = std::min(y0 + 1, (float)fNumY - 2);

			auto sx = 1.0 - tx;
			auto sy = 1.0 - ty;

			if (x0 < fNumX && y0 < fNumY)
				particleDensity[x0 * n + y0] += sx * sy;
			if (x1 < fNumX && y0 < fNumY)
				particleDensity[x1 * n + y0] += tx * sy;
			if (x1 < fNumX && y1 < fNumY)
				particleDensity[x1 * n + y1] += tx * ty;
			if (x0 < fNumX && y1 < fNumY)
				particleDensity[x0 * n + y1] += sx * ty;
		}

		if (particleRestDensity == 0.0) {
			auto sum           = 0.0;
			auto numFluidCells = 0;

			for (auto i = 0; i < fNumCells; i++) {
				if (cellType[i] == constants::FLUID_CELL) {
					sum += particleDensity[i];
					numFluidCells++;
				}
			}

			if (numFluidCells > 0)
				particleRestDensity = sum / numFluidCells;
		}

		// 			for (auto xi = 1; xi < this.fNumX; xi++) {
		// 				for (auto yi = 1; yi < this.fNumY; yi++) {
		// 					auto cellNr = xi * n + yi;
		// 					if (this.cellType[cellNr] != FLUID_CELL)
		// 						continue;
		// 					auto hx = this.h;
		// 					auto hy = this.h;

		// 					if (this.cellType[(xi - 1) * n + yi] == SOLID_CELL || this.cellType[(xi + 1) * n + yi] == SOLID_CELL)
		// 						hx -= this.particleRadius;
		// 					if (this.cellType[xi * n + yi - 1] == SOLID_CELL || this.cellType[xi * n + yi + 1] == SOLID_CELL)
		// 						hy -= this.particleRadius;

		// 					auto scale = this.h * this.h / (hx * hy)
		// 					d[cellNr] *= scale;
		// 				}
		// 			}
	}

	void transferVelocities(bool toGrid, float flipRatio)
	{
		auto n   = fNumY;
		auto h1  = fInvSpacing;
		float h2 = 0.5 * h;

		if (toGrid) {
			prevU = u;
			prevV = v;

			std::fill(du.begin(), du.end(), 0.0f);
			std::fill(dv.begin(), dv.end(), 0.0f);
			std::fill(u.begin(), u.end(), 0.0f);
			std::fill(v.begin(), v.end(), 0.0f);

			for (auto i = 0; i < fNumCells; i++)
				cellType[i] = s[i] == 0.0 ? constants::SOLID_CELL : constants::AIR_CELL;

			for (auto i = 0; i < numParticles; i++) {
				auto x      = particlePos[2 * i];
				auto y      = particlePos[2 * i + 1];
				auto xi     = std::clamp(floor(x * h1), 0.0f, (float)fNumX - 1);
				auto yi     = std::clamp(floor(y * h1), 0.0f, (float)fNumY - 1);
				auto cellNr = xi * n + yi;
				if (cellType[cellNr] == constants::AIR_CELL)
					cellType[cellNr] = constants::FLUID_CELL;
			}
		}

		for (auto component = 0; component < 2; component++) {
			auto dx = component == 0 ? 0.0f : h2;
			auto dy = component == 0 ? h2 : 0.0f;

			std::vector<float> f     = component == 0 ? u : v;
			std::vector<float> prevF = component == 0 ? prevU : prevV;
			std::vector<float> d     = component == 0 ? du : dv;

			for (auto i = 0; i < numParticles; i++) {
				auto x = particlePos[2 * i];
				auto y = particlePos[2 * i + 1];

				x = std::clamp(x, h, (fNumX - 1) * h);
				y = std::clamp(y, h, (fNumY - 1) * h);

				auto x0 = std::min(floor((x - dx) * h1), (float)fNumX - 2);
				auto tx = ((x - dx) - x0 * h) * h1;
				auto x1 = std::min(x0 + 1, (float)fNumX - 2);

				auto y0 = std::min(floor((y - dy) * h1), (float)fNumY - 2);
				auto ty = ((y - dy) - y0 * h) * h1;
				auto y1 = std::min(y0 + 1, (float)fNumY - 2);

				auto sx = 1.0 - tx;
				auto sy = 1.0 - ty;

				auto d0 = sx * sy;
				auto d1 = tx * sy;
				auto d2 = tx * ty;
				auto d3 = sx * ty;

				auto nr0 = x0 * n + y0;
				auto nr1 = x1 * n + y0;
				auto nr2 = x1 * n + y1;
				auto nr3 = x0 * n + y1;

				if (toGrid) {
					auto pv = particleVel[2 * i + component];
					f[nr0] += pv * d0;
					d[nr0] += d0;
					f[nr1] += pv * d1;
					d[nr1] += d1;
					f[nr2] += pv * d2;
					d[nr2] += d2;
					f[nr3] += pv * d3;
					d[nr3] += d3;
				} else {
					auto offset = component == 0 ? n : 1;
					auto valid0 = cellType[nr0] != constants::AIR_CELL || cellType[nr0 - offset] != constants::AIR_CELL ? 1.0 : 0.0;
					auto valid1 = cellType[nr1] != constants::AIR_CELL || cellType[nr1 - offset] != constants::AIR_CELL ? 1.0 : 0.0;
					auto valid2 = cellType[nr2] != constants::AIR_CELL || cellType[nr2 - offset] != constants::AIR_CELL ? 1.0 : 0.0;
					auto valid3 = cellType[nr3] != constants::AIR_CELL || cellType[nr3 - offset] != constants::AIR_CELL ? 1.0 : 0.0;

					auto v = particleVel[2 * i + component];
					auto d = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

					if (d > 0.0) {
						auto picV = (valid0 * d0 * f[nr0] + valid1 * d1 * f[nr1] + valid2 * d2 * f[nr2] + valid3 * d3 * f[nr3]) / d;
						auto corr = (valid0 * d0 * (f[nr0] - prevF[nr0]) + valid1 * d1 * (f[nr1] - prevF[nr1]) + valid2 * d2 * (f[nr2] - prevF[nr2]) +
						                valid3 * d3 * (f[nr3] - prevF[nr3])) /
						            d;
						auto flipV = v + corr;

						particleVel[2 * i + component] = (1.0 - flipRatio) * picV + flipRatio * flipV;
					}
				}
			}

			if (toGrid) {
				for (auto i = 0; i < f.size(); i++) {
					if (d[i] > 0.0)
						f[i] /= d[i];
				}

				// restore solid cells

				for (auto i = 0; i < fNumX; i++) {
					for (auto j = 0; j < fNumY; j++) {
						auto solid = cellType[i * n + j] == constants::SOLID_CELL;
						if (solid || (i > 0 && cellType[(i - 1) * n + j] == constants::SOLID_CELL))
							u[i * n + j] = prevU[i * n + j];
						if (solid || (j > 0 && cellType[i * n + j - 1] == constants::SOLID_CELL))
							v[i * n + j] = prevV[i * n + j];
					}
				}
			}
		}
	}

	void solveIncompressibility(int numIters, float dt, float overRelaxation, bool compensateDrift = true)
	{
		std::fill(p.begin(), p.end(), 0.0f);
		prevU = u;
		prevV = v;

		auto n  = fNumY;
		auto cp = density * h / dt;

		//for (auto i = 0; i < fNumCells; i++) {
		//	auto u = u[i]; ///eee pkoi^
		//	auto v = v[i];
		//}

		for (auto iter = 0; iter < numIters; iter++) {
			for (auto i = 1; i < fNumX - 1; i++) {
				for (auto j = 1; j < fNumY - 1; j++) {
					if (cellType[i * n + j] != constants::FLUID_CELL)
						continue;

					auto center = i * n + j;
					auto left   = (i - 1) * n + j;
					auto right  = (i + 1) * n + j;
					auto bottom = i * n + j - 1;
					auto top    = i * n + j + 1;

					//auto s   = s[center];
					auto sx0 = s[left];
					auto sx1 = s[right];
					auto sy0 = s[bottom];
					auto sy1 = s[top];
					auto s   = sx0 + sx1 + sy0 + sy1;
					if (s == 0.0)
						continue;

					auto div = u[right] - u[center] + v[top] - v[center];

					if (particleRestDensity > 0.0 && compensateDrift) {
						auto k           = 1.0;
						auto compression = particleDensity[i * n + j] - particleRestDensity;
						if (compression > 0.0)
							div = div - k * compression;
					}

					auto pValue = -div / s;
					pValue *= overRelaxation;
					p[center] += cp * pValue;

					u[center] -= sx0 * pValue;
					u[right] += sx1 * pValue;
					v[center] -= sy0 * pValue;
					v[top] += sy1 * pValue;
				}
			}
		}
	}

	void updateParticleColors()
	{
		// for (auto i = 0; i < this.numParticles; i++) {
		// 	this.particleColor[3 * i] *= 0.99;
		// 	this.particleColor[3 * i + 1] *= 0.99
		// 	this.particleColor[3 * i + 2] =
		// 		clamp(this.particleColor[3 * i + 2] + 0.001, 0.0, 1.0)
		// }

		// return;

		auto h1 = fInvSpacing;

		for (auto i = 0; i < numParticles; i++) {
			auto s = 0.01f;

			particleColor[3 * i]     = std::clamp(particleColor[3 * i] - s, 0.0f, 1.0f);
			particleColor[3 * i + 1] = std::clamp(particleColor[3 * i + 1] - s, 0.0f, 1.0f);
			particleColor[3 * i + 2] = std::clamp(particleColor[3 * i + 2] + s, 0.0f, 1.0f);

			auto x      = particlePos[2 * i];
			auto y      = particlePos[2 * i + 1];
			auto xi     = std::clamp(floor(x * h1), 1.0f, (float)fNumX - 1);
			auto yi     = std::clamp(floor(y * h1), 1.0f, (float)fNumY - 1);
			auto cellNr = xi * fNumY + yi;

			auto d0 = particleRestDensity;

			if (d0 > 0.0) {
				auto relDensity = particleDensity[cellNr] / d0;
				if (relDensity < 0.7) {
					auto s                   = 0.8f;
					particleColor[3 * i]     = s;
					particleColor[3 * i + 1] = s;
					particleColor[3 * i + 2] = 1.0f;
				}
			}
		}
	}

	void setSciColor(int cellNr, float val, float minVal, float maxVal)
	{
		val     = std::min(std::max(val, minVal), maxVal - 0.0001f);
		auto d  = maxVal - minVal;
		val     = d == 0.0f ? 0.5f : (val - minVal) / d;
		float m = 0.25f;
		int num = floor(val / m);
		float s = (val - num * m) / m;
		float r, g, b;

		switch (num) {
			case 0:
				r = 0.0f;
				g = s;
				b = 1.0f;
				break;
			case 1:
				r = 0.0f;
				g = 1.0f;
				b = 1.0f - s;
				break;
			case 2:
				r = s;
				g = 1.0f;
				b = 0.0f;
				break;
			case 3:
				r = 1.0f;
				g = 1.0f - s;
				b = 0.0f;
				break;
		}

		cellColor[3 * cellNr]     = r;
		cellColor[3 * cellNr + 1] = g;
		cellColor[3 * cellNr + 2] = b;
	}

	void updateCellColors()
	{
		std::fill(cellColor.begin(), cellColor.end(), 0.5f);

		for (auto i = 0; i < fNumCells; i++) {
			if (cellType[i] == constants::SOLID_CELL) {
				cellColor[3 * i]     = 0.5f;
				cellColor[3 * i + 1] = 0.5f;
				cellColor[3 * i + 2] = 0.5f;
			} else if (cellType[i] == constants::FLUID_CELL) {
				auto d = particleDensity[i];
				if (particleRestDensity > 0.0f)
					d /= particleRestDensity;
				setSciColor(i, d, 0.0f, 2.0f);
			}
		}
	}

	void simulate(float dt,
	    float gravity,
	    float flipRatio,
	    int numPressureIters,
	    int numParticleIters,
	    float overRelaxation,
	    bool compensateDrift,
	    bool separateParticles,
	    int obstacleX,
	    int obstacleY,
	    float obstacleRadius)
	{
		auto numSubSteps = 1;
		auto sdt         = dt / numSubSteps;

		for (auto step = 0; step < numSubSteps; step++) {
			integrateParticles(sdt, gravity);
			if (separateParticles)
				pushParticlesApart(numParticleIters);
			handleParticleCollisions(obstacleX, obstacleY, obstacleRadius);
			transferVelocities(true, flipRatio);
			updateParticleDensity();
			solveIncompressibility(numPressureIters, sdt, overRelaxation, compensateDrift);
			transferVelocities(true, flipRatio);
		}

		updateParticleColors();
		updateCellColors();
	}
};



namespace Ui {
class fluid;
}

class fluid : public QDialog
{
	Q_OBJECT

public:
	explicit fluid(QWidget* parent = nullptr);

	~fluid() { delete ui; }


	void setupScene();
	void update();


protected:
private slots:
	void on_pauseBtn_clicked()
	{
		m_paused = !m_paused;
		update();
	}
	void on_resetBtn_clicked() { setupScene(); }

	void on_dsBox4_valueChanged(double d)
	{
		constants::scale    = d;
		constants::simwidth = (float)constants::maxwidth / constants::scale;
	}


private:
	Ui::fluid* ui;
	QGraphicsScene* m_scene;

	QImage m_imageCanvas;

	QGraphicsPixmapItem* m_sceneItem;
	std::vector<QGraphicsRectItem*> m_gridItems;
	std::vector<QGraphicsEllipseItem*> m_particleItems;
	QGraphicsEllipseItem* m_obstacleItem;

	int m_frameNr;
	float m_obstacleX;
	float m_obstacleY;
	bool m_paused;

	FlipFluid m_f;


	void setupObstacle(int x, int y, bool reset);

	void simulate();
	void draw();
	void requestAnimationFrame();



	QElapsedTimer m_timer;
};

#endif