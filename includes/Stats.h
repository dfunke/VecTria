#pragma once

#include <unordered_map>

class Stats {

public:
    static Stats &getInstance() {
        static Stats instance;

        return instance;
    }

private:
    Stats() {}

public:
    Stats(const Stats &) = delete;

    void operator=(const Stats &) = delete;

private:
    std::unordered_map<std::string, std::pair<double, std::size_t>> m_values;

public:
    void clear() {
        m_values.clear();
    }

    void incStat(const std::string &key, const double val = 1) {
        auto it = m_values.find(key);
        if (it == m_values.end()) {
            m_values[key] = std::make_pair(0, 0);
        }

        m_values[key].first += val;
        m_values[key].second += 1;
    }

    void updateStat(const std::string &key, const double val = 1) {
        auto it = m_values.find(key);
        if (it == m_values.end()) {
            m_values[key] = std::make_pair(0, 0);
        }

        m_values[key].first = m_values[key].first + (val - m_values[key].first) / (m_values[key].second + 1);
        m_values[key].second += 1;
    }

    double getStat(const std::string &key) const {
        auto it = m_values.find(key);
        if (it == m_values.end()) {
            return 0;
        }

        return it->second.first;
    }

};

#ifdef ENABLE_STATS

#define STAT_INC(stat) Stats::getInstance().incStat(#stat)
#define STAT_ADD(stat, val) Stats::getInstance().incStat(#stat, val)
#define STAT_UPD(stat, val) Stats::getInstance().updateStat(#stat, val)
#define STAT_GET(stat) Stats::getInstance().getStat(#stat)
#define STAT_CLEAR Stats::getInstance().clear()

#else // ENABLE_STATS

#define STAT_INC(stat) ((void)(0))
#define STAT_GET 0
#define STAT_ADD(stat, val) ((void)(0))
#define STAT_UPD(stat, val) ((void)(0))
#define STAT_CLEAR ((void)(0))

#endif