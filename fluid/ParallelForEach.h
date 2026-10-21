#pragma once

#include <algorithm>
#include <iterator>
#include <utility>
#include <QtConcurrent/QtConcurrentMap>

template <typename Iterator, typename Function>
void fluidParallelForEach(Iterator first, Iterator last, Function function)
{
	// Avoid scheduling overhead for small ranges. QtConcurrent uses the
	// configured Qt thread pool on both native and threaded-WASM builds.
	if (std::distance(first, last) < 1024) {
		std::for_each(first, last, function);
		return;
	}

	QtConcurrent::blockingMap(first, last, std::move(function));
}
