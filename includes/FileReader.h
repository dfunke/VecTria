#pragma once

#include <fstream>
#include <filesystem>
#include <vector>

#include <boost/algorithm/string.hpp>

#include "GeometryStructures.h"

template<tDimType D, typename Precision>
struct FileReader;

template<typename Precision>
struct FileReader<3, Precision> {

    template<class PointArray>
    tIndexType read(const std::string &filename, PointArray &pa) const {

        std::ifstream f(filename);
        assert(f.is_open());
        std::string line;

        // first item is the timestep
        std::getline(f, line);
        assert(line == "ITEM: TIMESTEP");

        std::getline(f, line);
        tIndexType timestep = static_cast<tIndexType>(std::stol(line));

        // the timestep data is weird, use filename instead for now
        timestep = static_cast<tIndexType>(stol(std::filesystem::path(filename).stem().stem().string()));

        // second item is the number of atoms
        std::getline(f, line);
        assert (line == "ITEM: NUMBER OF ATOMS");
        std::getline(f, line);
        tIndexType n = static_cast<tIndexType>(stol(line));


        // third item is the bounding box
        std::getline(f, line);
        assert (line.find("ITEM: BOX BOUNDS") != std::string::npos);

        std::getline(f, line);
        //bbox[0, :] = np.asarray(f.readline().split(' '), dtype=np.float)
        std::getline(f, line);
        //bbox[1, :] = np.asarray(f.readline().split(' '), dtype=np.float)
        std::getline(f, line);
        //bbox[2, :] = np.asarray(f.readline().split(' '), dtype=np.float)


        // fourth item are the atoms
        std::getline(f, line);
        assert (line.find("ITEM: ATOMS") != std::string::npos);
        pa.ensure(n);

        Precision maxx = 0, maxy = 0, maxz = 0;
        for (tIndexType i = 0; i < n; ++i) {
            std::getline(f, line);

            std::vector<std::string> splits;
            boost::split(splits, line, [](char c) { return c == ' '; });

            auto p = pa.get(static_cast<tIndexType>(std::stol(splits[0]) - 1));

            p[0] = static_cast<Precision>(std::stod(splits[2]));
            p[1] = static_cast<Precision>(std::stod(splits[3]));
            p[2] = static_cast<Precision>(std::stod(splits[4]));

            if(p[0] > maxx) maxx = p[0];
            if(p[1] > maxy) maxy = p[1];
            if(p[2] > maxz) maxz = p[2];
        }

        Predicates<Precision>::set_static_limits(maxx, maxy, maxz);

        return timestep;
    }

    template<class PointArray>
    PointArray read(const std::string &filename) const {
        PointArray pa;
        read(filename, pa);
        return pa;
    }
};
