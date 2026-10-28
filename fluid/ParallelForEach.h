#pragma once

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <numeric>
#include <utility>
#include <vector>
#include <QtConcurrent/QtConcurrentMap>

template <typename Iterator, typename Function>
void fluidParallelForEach(Iterator first, Iterator last, Function function)
{
	// Avoid scheduling overhead for small ranges. QtConcurrent uses the
	// configured Qt thread pool on both native and threaded-WASM builds.
	if (std::distance(first, last) < 512) {
		std::for_each(first, last, function);
		return;
	}

	QtConcurrent::blockingMap(first, last, std::move(function));
}

template <typename Function>
void fluidParallelForEachIndex(const std::size_t count, Function function)
{
	if (count < 512) {
		for (std::size_t index = 0; index < count; ++index)
			function(index);
		return;
	}

	// Map blocks rather than individual indices. The inner loop gives the
	// optimizer a contiguous range over the particle SoA arrays to vectorize.
	constexpr std::size_t blockSize = 512;
	const std::size_t blockCount = (count + blockSize - 1) / blockSize;
	std::vector<std::size_t> blocks(blockCount);
	std::iota(blocks.begin(), blocks.end(), std::size_t { 0 });
	QtConcurrent::blockingMap(blocks.begin(), blocks.end(), [&](const std::size_t block) {
		const std::size_t first = block * blockSize;
		const std::size_t last = std::min(first + blockSize, count);
		for (std::size_t index = first; index < last; ++index)
			function(index);
	});
}
