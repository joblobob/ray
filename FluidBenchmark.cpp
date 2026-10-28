#include "FluidBenchmark.h"

#include <QCoreApplication>
#include <QDebug>
#include <QDir>
#include <QElapsedTimer>
#include <QFile>
#include <QFileInfo>
#include <QTextStream>

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <string>
#include <vector>

#include <mdspan>

import BaseStructures;
import Constants;
import FlipFluid;
import Obstacle;

namespace {
struct Samples {
	std::vector<double> microseconds;

	void add(const double ns) { microseconds.push_back(ns * 0.001); }
};

void initializeBenchmarkScene(FlipFluid& fluid)
{
	ParticleStore& particles = fluid.particleMap;
	for (int i = 0; i < constants::numX; ++i) {
		for (int j = 0; j < constants::numY; ++j) {
			const std::size_t id = static_cast<std::size_t>(i * constants::numY + j);
			particles.posX[id] = constants::particleStartX + constants::particleSpacingX * i + (j % 2 == 0 ? 0.0f : constants::particleStaggerX);
			particles.posY[id] = constants::particleStartY + constants::particleRowSpacing * j;
			particles.velX[id] = 0.0f;
			particles.velY[id] = 0.0f;
			particles.colorR[id] = 0.0f;
			particles.colorG[id] = 0.0f;
			particles.colorB[id] = 1.0f;
		}
	}

	std::mdspan grid(fluid.gridCells.data(), constants::fNumX, constants::fNumY);
	for (std::size_t i = 0; i < grid.extent(0); ++i) {
		for (std::size_t j = 0; j < grid.extent(1); ++j) {
			Cell& cell = grid[i, j];
			cell.s = (i == 0 || i + 1 == grid.extent(0) || j == 0) ? 0.0f : 1.0f;
			cell.cellNumX = static_cast<int>(i);
			cell.cellNumY = static_cast<int>(j);
		}
	}

	// A fixed, submerged obstacle gives each run the same interaction workload.
	setupObstacle(fluid.gridCells, 300.0f, 200.0f, fluid.obstacle, true);
}

bool readPositiveInt(const QStringList& args, const QString& option, int& value)
{
	const qsizetype index = args.indexOf(option);
	if (index < 0)
		return true;
	if (index + 1 >= args.size())
		return false;
	bool ok = false;
	const int parsed = args[index + 1].toInt(&ok);
	if (!ok || parsed < 1)
		return false;
	value = parsed;
	return true;
}

} // namespace

int runFluidBenchmark(const QStringList& arguments)
{
	int warmupFrames = 30;
	int measuredFrames = 180;
	if (!readPositiveInt(arguments, QStringLiteral("--warmup"), warmupFrames)
	    || !readPositiveInt(arguments, QStringLiteral("--frames"), measuredFrames)) {
		qCritical().noquote() << "Usage: ray1 --benchmark [--warmup N] [--frames N] [--output FILE]";
		return 2;
	}

	QString outputPath = QCoreApplication::applicationDirPath() + QStringLiteral("/fluid_benchmark.csv");
	const qsizetype outputOption = arguments.indexOf(QStringLiteral("--output"));
	if (outputOption >= 0) {
		if (outputOption + 1 >= arguments.size()) {
			qCritical().noquote() << "--output requires a file path";
			return 2;
		}
		const QFileInfo requested(arguments[outputOption + 1]);
		outputPath = requested.isAbsolute() ? requested.absoluteFilePath() : QDir::current().absoluteFilePath(requested.filePath());
	}

	FlipFluid fluid(constants::maxwidth, constants::maxheight, constants::cellHeight, constants::maxParticles);
	initializeBenchmarkScene(fluid);

	for (int i = 0; i < warmupFrames; ++i)
		fluid.simulate(constants::dt, false);

	std::map<std::string, Samples> samples;
	for (int frame = 0; frame < measuredFrames; ++frame) {
		QElapsedTimer frameTimer;
		frameTimer.start();
		const std::vector<ExecutionLog> stages = fluid.simulate(constants::dt, true);
		samples["simulate_total"].add(static_cast<double>(frameTimer.nsecsElapsed()));
		for (const ExecutionLog& stage : stages)
			samples[stage.message].add(static_cast<double>(stage.nsElapsed));
	}

	QFile csvFile(outputPath);
	if (!csvFile.open(QIODevice::WriteOnly | QIODevice::Truncate | QIODevice::Text)) {
		qCritical().noquote() << "Could not write benchmark CSV:" << outputPath;
		return 1;
	}
	QTextStream csv(&csvFile);
	csv.setRealNumberPrecision(9);
	csv << "metric,samples,mean_us,median_us,p95_us,min_us,max_us,particles,grid_x,grid_y,dt_seconds,warmup_frames\n";
	for (auto& [name, series] : samples) {
		auto& values = series.microseconds;
		std::sort(values.begin(), values.end());
		const double mean = std::accumulate(values.begin(), values.end(), 0.0) / static_cast<double>(values.size());
		const double median = values[values.size() / 2];
		const std::size_t p95Index = static_cast<std::size_t>(std::ceil(values.size() * 0.95)) - 1;
		csv << QString::fromStdString(name) << ',' << values.size() << ',' << mean << ',' << median << ',' << values[p95Index] << ',' << values.front() << ',' << values.back() << ','
		    << constants::maxParticles << ',' << constants::fNumX << ',' << constants::fNumY << ',' << constants::dt << ',' << warmupFrames << '\n';
	}
	csv.flush();
	csvFile.close();

	qInfo().noquote() << "Benchmark complete:" << outputPath;
	qInfo().noquote() << "particles=" << constants::maxParticles << " grid=" << constants::fNumX << 'x' << constants::fNumY
	                  << " warmup_frames=" << warmupFrames << " measured_frames=" << measuredFrames;
	return 0;
}
