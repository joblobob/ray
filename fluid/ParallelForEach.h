#pragma once

#include <QtConcurrent/QtConcurrentMap>
#include <algorithm>
#include <cstddef>
#include <iterator>
#include <numeric>
#include <utility>
#include <vector>

template<typename Iterator, typename Function> void fluidParallelForEach(Iterator first, Iterator last, Function function)
{
	// Avoid scheduling overhead for small ranges. QtConcurrent uses the
	// configured Qt thread pool on both native and threaded-WASM builds.
	auto n = std::distance(first, last);
	if (n < 512) {
		std::for_each(first, last, function);
		return;
	}

	QtConcurrent::blockingMap(first, last, std::move(function));
}

template<typename Function> void fluidParallelForEachIndex(const std::size_t count, Function function)
{
	if (count < 512) {
		for (std::size_t index = 0; index < count; ++index)
			function(index);
		return;
	}

	// Map blocks rather than individual indices. Larger block size reduces
	// QtConcurrent task overhead and improves cache reuse.
	constexpr std::size_t blockSize = 1024;
	const std::size_t blockCount    = (count + blockSize - 1) / blockSize;

	// Reuse a thread_local buffer to avoid per-call allocation.
	thread_local std::vector<std::size_t> blocks;
	blocks.resize(blockCount);
	std::iota(blocks.begin(), blocks.end(), std::size_t { 0 });

	QtConcurrent::blockingMap(blocks.begin(), blocks.end(), [&](const std::size_t block) {
		const std::size_t first = block * blockSize;
		const std::size_t last  = std::min(first + blockSize, count);
		for (std::size_t index = first; index < last; ++index)
			function(index);
	});
}
