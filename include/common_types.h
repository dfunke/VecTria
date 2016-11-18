#pragma once

using tIdType = unsigned int;

using tIdSet = LP_Set<tIdType, true>;
using tConcurrentFixedIdSet<tIdType, false>;
using tConcurrentGrowingIdSet = GrowingHashTable <Concurrent_LP_Set<tIdType, true>>;
using tConcurrentGrowingIdSetHandle = GrowingHashTableHandle <Concurrent_LP_Set<tIdType, true>>;