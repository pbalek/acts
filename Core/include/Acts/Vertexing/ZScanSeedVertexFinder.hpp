// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"

#include <vector>
#include <string>

namespace Acts {
/// @class ZScanSeedVertexFinder
///
/// @brief Implements the vertex finder based on the track seeds

template <typename spacepoint_t>
class ZScanSeedVertexFinder {
    public:
        /// Configuration struct
        struct Config{
            /// @brief Config constructor
            ///
            /// @param maxZRdev maximum deviation in z-r plane
            /// @param maxXYdev maximum deviation in x-y plane
            // Config() {}
            // Config(int maxZRdev, int maxXYdev) : maxZRdeviation(maxZRdev), maxXYdeviation(maxXYdev) {}

            // maximum deviation in z-r plane between the first 2 spacepoints and the last 2 spacepoints
            float maxZRdeviation = 0.1;
            // maximum deviation in z-r plane between the first 2 spacepoints and the last 2 spacepoints
            float maxXYdeviation = 0.1;

            // maximum distance in r between two measurements within one seed
            // float deltaRMax = 270 * Acts::UnitConstants::mm;
            // range defined by rMinMiddle and rMaxMiddle
            float rMinMiddle = 60.f * Acts::UnitConstants::mm;
            float rMaxMiddle = 120.f * Acts::UnitConstants::mm;

            float zBinSize = 1.f * Acts::UnitConstants::mm;

            /// We allow multiple space point collections to allow different parts of
            /// the detector to use different algorithms for space point construction,
            /// e.g. single-hit space points for pixel-like detectors or double-hit
            /// space points for strip-like detectors.
            // std::vector<std::string> inputSpacePoints;
            /// Output vertex collection.
            // std::string outputVertices;
        };
        
        /// Const access to the config
        const Config& config() const { return m_cfg; }

        ZScanSeedVertexFinder(Config& cfg);
        ZScanSeedVertexFinder() {Config cfg; m_cfg = cfg;}
        ~ZScanSeedVertexFinder() = default;

        // TODO: delete default constructor etc.?

        std::vector<float> findVertex(const std::vector<spacepoint_t>& spacepoints) const;

    private:
        struct Triplet{
            Triplet(spacepoint_t aa, spacepoint_t bb, spacepoint_t cc) : a(aa), b(bb), c(cc) {}
            spacepoint_t a,b,c;
        };

        Acts::ZScanSeedVertexFinder<spacepoint_t>::Config m_cfg;

        std::vector<std::vector<spacepoint_t>> sortSpacepoints(const std::vector<spacepoint_t>& spacepoints) const;

        std::vector<typename Acts::ZScanSeedVertexFinder<spacepoint_t>::Triplet> findTriplets(const std::vector<std::vector<spacepoint_t>>& sorted_spacepoints) const;
        bool isTripletValid(const Acts::ZScanSeedVertexFinder<spacepoint_t>::triplet) const;

        std::vector<int> makeZHist(const std::vector<Acts::ZScanSeedVertexFinder<spacepoint_t>::Triplet>& triplets) const;
        
        float findZPeak(const std::vector<int>& hist) const;
        
};

} // namespace Acts

#include "Acts/Vertexing/ZScanSeedVertexFinder.ipp"
