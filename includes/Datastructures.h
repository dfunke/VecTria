#pragma once

#include <vector>
#include <array>

using tIndexType = unsigned long;
constexpr tIndexType INF = ~tIndexType(0);

using tDimType = unsigned short;

constexpr tDimType X = 0;
constexpr tDimType Y = 1;
constexpr tDimType Z = 2;

template<typename Precision>
using tFloatVector = std::vector<Precision>;
using tIndexVector = std::vector<tIndexType>;

template<typename T, tDimType D>
class MemoryLayoutPA : private std::vector<T> {
private:
    using base = std::vector<T>;

public:

    inline const T &operator()(const tIndexType &i, const tDimType &d) const {
        return base::operator[](D * i + d);
    }

    inline T &operator()(const tIndexType &i, const tDimType &d) {
        return base::operator[](D * i + d);
    }

    inline auto size() const {
        return base::size() / D;
    }

    inline void ensure(const tIndexType &i) {
        if (base::size() < D * (i + 1)) {
            base::resize(D * (i + 1));
        }
    }

    inline void ensure(const tIndexType &i) const {
        if (base::size() < D * (i + 1)) {
            throw std::out_of_range(
                    "index " + std::to_string(i) + " is out of range " + std::to_string(size() / D));
        }
    }

};

template<typename T, tDimType D>
class MemoryLayoutAoA : private std::array<std::vector<T>, D> {
private:
    using base = std::array<std::vector<T>, D>;

public:

    inline const T &operator()(const tIndexType &i, const tDimType &d) const {
        return base::operator[](d)[i];
    }

    inline T &operator()(const tIndexType &i, const tDimType &d) {
        return base::operator[](d)[i];
    }

    inline auto size() const {
        return base::operator[](0).size();
    }

    inline void ensure(const tIndexType &i) {
        for (tDimType d = 0; d < D; ++d) {
            if (base::operator[](d).size() < i + 1) {
                base::operator[](d).resize(i + 1);
            }
        }
    }

    inline void ensure(const tIndexType &i) const {
        for (tDimType d = 0; d < D; ++d) {
            if (base::operator[](d).size() < i + 1) {
                throw std::out_of_range(
                        "index " + std::to_string(i) + " is out of range " + std::to_string(size()));
            }
        }
    }

};