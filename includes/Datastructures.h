#pragma once

#include <vector>
#include <array>

#ifdef HAS_Vc

#include <Vc/Vc>

#endif

using tIndexType = unsigned int;
constexpr tIndexType INF = ~tIndexType(0);

using tDimType = unsigned short;

constexpr tDimType X = 0;
constexpr tDimType Y = 1;
constexpr tDimType Z = 2;

#ifndef NDEBUG
#define ACCESS at
#else
#define ACCESS operator[]
#endif

template<typename Precision>
using tFloatVector = std::vector<Precision>;
using tIndexVector = std::vector<tIndexType>;

template<typename T, tDimType D>
class MemoryLayoutPA : private std::vector<T> {
private:
    using base = std::vector<T>;
    tIndexType m_size = 0;

public:
    static constexpr bool isVectorized = false;

public:

    inline const T &operator()(const tIndexType &i, const tDimType &d) const {
        return base::ACCESS(D * i + d);
    }

    inline T &operator()(const tIndexType &i, const tDimType &d) {
        return base::ACCESS(D * i + d);
    }

    inline auto size() const {
        return m_size;
    }

    inline void ensure(const tIndexType &i) {
        if (m_size < i) {
            base::resize(D * (i + 1));
            m_size = i + 1;
        }
    }

    inline void ensure(const tIndexType &i) const {
        if (m_size < i) {
            throw std::out_of_range(
                    "index " + std::to_string(i) + " is out of range " + std::to_string(size() / D));
        }
    }

    static std::string name() {
        return "PA";
    }

};

template<typename T, tDimType D>
class MemoryLayoutAoA : private std::array<std::vector<T>, D> {
private:
    using base = std::array<std::vector<T>, D>;
    tIndexType m_size = 0;

public:
    static constexpr bool isVectorized = false;

public:

    inline const T &operator()(const tIndexType &i, const tDimType &d) const {
        return base::ACCESS(d).ACCESS(i);
    }

    inline T &operator()(const tIndexType &i, const tDimType &d) {
        return base::ACCESS(d).ACCESS(i);
    }

    inline auto size() const {
        return m_size;
    }

    inline void ensure(const tIndexType &i) {
        if (m_size < i) {
            for (tDimType d = 0; d < D; ++d) {
                base::ACCESS(d).resize(i + 1);
            }
            m_size = i + 1;
        }
    }

    inline void ensure(const tIndexType &i) const {
        if (m_size < i) {
            throw std::out_of_range(
                    "index " + std::to_string(i) + " is out of range " + std::to_string(size()));
        }
    }

    static std::string name() {
        return "AoA";
    }

};

#ifdef HAS_Vc

template<typename T, tDimType D>
class MemoryLayoutVectorizedPA : private std::vector<T> {
private:
    using base = std::vector<T>;
    tIndexType m_size = 0;

public:
    static constexpr bool isVectorized = true;

public:

    inline const T &operator()(const tIndexType &i, const tDimType &d) const {
        return base::ACCESS(D * i + d);
    }

    inline T &operator()(const tIndexType &i, const tDimType &d) {
        return base::ACCESS(D * i + d);
    }

    inline Vc::Vector<T> vec(const Vc::Vector<tIndexType> &i, const tDimType &d) const {
        static_assert(Vc::Vector<T>::size() == Vc::Vector<tIndexType>::size());
        return Vc::Vector<T>(base::data(), Vc::Vector<tIndexType>(D) * i + Vc::Vector<tIndexType>(d));
    }

    inline Vc::Vector<T> vec(const tIndexType &i, const tDimType &d) const {
        static_assert(Vc::Vector<T>::size() == Vc::Vector<tIndexType>::size());
        return Vc::Vector<T>(base::data(), Vc::Vector<tIndexType>(D) *
                                           (Vc::Vector<tIndexType>(i) + Vc::Vector<tIndexType>::IndexesFromZero()) +
                                           Vc::Vector<tIndexType>(d));
    }

    inline Vc::Vector<T> vec(const tIndexType &i) const {
        static_assert(D <= Vc::Vector<T>::size());
        return Vc::Vector<T>(base::data()[D * i]);
    }

    inline void store(const Vc::Vector<tIndexType> &i, const tDimType &d, const Vc::Vector<T> &data) {
        static_assert(Vc::Vector<T>::size() == Vc::Vector<tIndexType>::size());
        data.scatter(base::data(), Vc::Vector<tIndexType>(D) * i + Vc::Vector<tIndexType>(d));
    }

    inline void store(const tIndexType &i, const tDimType &d, const Vc::Vector<T> &data) {
        static_assert(Vc::Vector<T>::size() == Vc::Vector<tIndexType>::size());
        data.scatter(base::data(), Vc::Vector<tIndexType>(D) *
                                   (Vc::Vector<tIndexType>(i) + Vc::Vector<tIndexType>::IndexesFromZero()) +
                                   Vc::Vector<tIndexType>(d));
    }

    inline auto size() const {
        return m_size;
    }

    inline void ensure(const tIndexType &i) {
        if (m_size < i) {
            base::resize(D * (i + 1) + (D < Vc::Vector<T>::size() ? Vc::Vector<T>::size() - D : 0));
            m_size = i + 1;
        }
    }

    inline void ensure(const tIndexType &i) const {
        if (m_size < i) {
            throw std::out_of_range(
                    "index " + std::to_string(i) + " is out of range " + std::to_string(size() / D));
        }
    }

    static std::string name() {
        return "vectorized PA";
    }

};

template<typename T, tDimType D>
class MemoryLayoutVectorizedGroupedPA : private std::vector<T> {
private:
    using base = std::vector<T>;
    tIndexType m_size = 0;

public:
    static constexpr bool isVectorized = true;

private:
    static constexpr tIndexType N = Vc::Vector<T>::size();
    static constexpr tIndexType ND = N * D;

public:

    inline const T &operator()(const tIndexType &i, const tDimType &d) const {
        return base::ACCESS((i / N) * ND + d * N + i % N);
    }

    inline T &operator()(const tIndexType &i, const tDimType &d) {
        return base::ACCESS((i / N) * ND + d * N + i % N);
    }

    inline Vc::Vector<T> vec(const Vc::Vector<tIndexType> &i, const tDimType &d) const {
        static_assert(Vc::Vector<T>::size() == Vc::Vector<tIndexType>::size());
        return Vc::Vector<T>(base::data(),
                             (i / Vc::Vector<tIndexType>(N)) * Vc::Vector<tIndexType>(ND)
                             + Vc::Vector<tIndexType>(d) * Vc::Vector<tIndexType>(N) +
                             i % Vc::Vector<tIndexType>(N));
    }

    inline Vc::Vector<T> vec(const tIndexType &i, const tDimType &d) const {
        static_assert(Vc::Vector<T>::size() == Vc::Vector<tIndexType>::size());
        return Vc::Vector<T>(&base::data()[(i / N) * ND + d * N]);
    }

    inline Vc::Vector<T> vec(const tIndexType &i) const {
        static_assert(D <= Vc::Vector<T>::size());
        return Vc::Vector<T>(base::data(), (i / Vc::Vector<tIndexType>(N)) * Vc::Vector<tIndexType>(ND)
                                           + Vc::Vector<tIndexType>::IndexesFromZero() * Vc::Vector<tIndexType>(N) +
                                           i % Vc::Vector<tIndexType>(N));
    }

    inline void store(const Vc::Vector<tIndexType> &i, const tDimType &d, const Vc::Vector<T> &data) {
        static_assert(Vc::Vector<T>::size() == Vc::Vector<tIndexType>::size());
        data.scatter(base::data(), (i / Vc::Vector<tIndexType>(N)) * Vc::Vector<tIndexType>(ND)
                                   + Vc::Vector<tIndexType>(d) * Vc::Vector<tIndexType>(N) +
                                   i % Vc::Vector<tIndexType>(N));
    }

    inline void store(const tIndexType &i, const tDimType &d, const Vc::Vector<T> &data) {
        static_assert(Vc::Vector<T>::size() == Vc::Vector<tIndexType>::size());
        data.store(&base::data()[(i / N) * ND + d * N]);
    }

    inline auto size() const {
        return m_size;
    }

    inline void ensure(const tIndexType &i) {
        if (m_size < i) {
            base::resize(((i / N) + 1) * ND + (D < Vc::Vector<T>::size() ? Vc::Vector<T>::size() - D : 0));
            m_size = i + 1;
        }
    }

    inline void ensure(const tIndexType &i) const {
        if (m_size < i) {
            throw std::out_of_range(
                    "index " + std::to_string(i) + " is out of range " + std::to_string(size() / D));
        }
    }

    static std::string name() {
        return "vectorized grouped PA";
    }

};

template<typename T, tDimType D>
class MemoryLayoutVectorizedAoA : private std::array<std::vector<T>, D> {
private:
    using base = std::array<std::vector<T>, D>;
    tIndexType m_size = 0;

public:
    static constexpr bool isVectorized = true;

public:

    inline const T &operator()(const tIndexType &i, const tDimType &d) const {
        return base::ACCESS(d).ACCESS(i);
    }

    inline T &operator()(const tIndexType &i, const tDimType &d) {
        return base::ACCESS(d).ACCESS(i);
    }

    inline Vc::Vector<T> vec(const Vc::Vector<tIndexType> &i, const tDimType &d) const {
        static_assert(Vc::Vector<T>::size() == Vc::Vector<tIndexType>::size());
        return Vc::Vector<T>(base::ACCESS(d).data(), i);
    }

    inline Vc::Vector<T> vec(const tIndexType &i, const tDimType &d) const {
        static_assert(Vc::Vector<T>::size() == Vc::Vector<tIndexType>::size());
        return Vc::Vector<T>(&base::ACCESS(d).data()[i]);
    }

    inline Vc::Vector<T> vec(const tIndexType &i) const {
        static_assert(D <= Vc::Vector<T>::size());
        Vc::Vector<T> v;
        for (tDimType d = 0; d < D; ++d) {
            v[d] = base::ACCESS(d).ACCESS(i);
        }

        return v;
    }

    inline void store(const Vc::Vector<tIndexType> &i, const tDimType &d, const Vc::Vector<T> &data) {
        static_assert(Vc::Vector<T>::size() == Vc::Vector<tIndexType>::size());
        data.scatter(base::ACCESS(d).data(), i);
    }

    inline void store(const tIndexType &i, const tDimType &d, const Vc::Vector<T> &data) {
        static_assert(Vc::Vector<T>::size() == Vc::Vector<tIndexType>::size());
        data.store(&base::ACCESS(d).data()[i]);
    }

    inline auto size() const {
        return m_size;
    }

    inline void ensure(const tIndexType &i) {
        if (m_size < i) {
            for (tDimType d = 0; d < D; ++d) {
                base::ACCESS(d).resize(i + 1);
            }
            m_size = i + 1;
        }
    }

    inline void ensure(const tIndexType &i) const {

        if (m_size < i) {
            throw std::out_of_range(
                    "index " + std::to_string(i) + " is out of range " + std::to_string(size()));
        }

    }

    static std::string name() {
        return "vectorized AoA";
    }

};

#endif

