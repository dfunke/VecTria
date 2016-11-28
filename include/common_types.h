#pragma once

#include <vector>
#include <cstdint>

#include <tbb/concurrent_vector.h>

using tIdType = uint64_t;
using tHashType = uint32_t;
using tMarkType = uint8_t;
using tDimType = uint8_t;
using tRecDepthType = uint8_t;

using tIdVector = std::vector<tIdType>;
using tConcurrentIdVector = tbb::concurrent_vector<tIdType>;

using tIdSet = LP_Set<tIdType, true>;
using tConcurrentFixedIdSet = Concurrent_LP_Set<tIdType, false>;
using tConcurrentGrowingIdSet = GrowingHashTable <Concurrent_LP_Set<tIdType, true>>;
using tConcurrentGrowingIdSetHandle = GrowingHashTableHandle <Concurrent_LP_Set<tIdType, true>>;

using tFacetteLookUp = LP_MultiMap<tHashType, tIdType>;
using tConcurrentGrowingFacetteLookUp = GrowingHashTable <Concurrent_LP_MultiMap<tHashType, tIdType>>;
using tConcurrentGrowingFacetteLookUpHandle= GrowingHashTableHandle <Concurrent_LP_MultiMap<tHashType, tIdType>>;