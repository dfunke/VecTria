#pragma once

#include <vector>

using tIndexType = unsigned long;
constexpr tIndexType INF = ~tIndexType(0);

using tDimType = unsigned short;

constexpr tDimType X = 0;
constexpr tDimType Y = 1;
constexpr tDimType Z = 2;

template<typename Precision>
using tFloatVector = std::vector<Precision>;
using tIndexVector = std::vector<tIndexType>;
