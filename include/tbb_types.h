#pragma once

#include <tbb/enumerable_thread_specific.h>

template<typename T>
using tbbETS = tbb::enumerable_thread_specific <T,
tbb::cache_aligned_allocator<T>,
tbb::ets_key_usage_type::ets_key_per_instance>;