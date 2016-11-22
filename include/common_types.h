#pragma once

using tIdType = uint64_t;
using tHashType = uint32_t;
using tMarkType = uint8_t;

using tIdSet = LP_Set<tIdType, true>;
using tConcurrentFixedIdSet<tIdType, false>;
using tConcurrentGrowingIdSet = GrowingHashTable <Concurrent_LP_Set<tIdType, true>>;
using tConcurrentGrowingIdSetHandle = GrowingHashTableHandle <Concurrent_LP_Set<tIdType, true>>;

using tConcurrentGrowingFacetteLookUp = GrowingHashTable<Concurrent_LP_MultiMap<tHashType, tIdType>>;
using tConcurrentGrowingFacetteLookUpHandle= GrowingHashTableHandle<Concurrent_LP_MultiMap<tHashType, tIdType>>;