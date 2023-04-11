// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <vector>
#include <string>

namespace Acts {
/// @class ZScanSeedVertexFinder
///
/// @brief Implements the vertex finder based on the track seeds
/// 0. Assumes there is only 1 vertex and that it has a high multiplicity
/// 1. Sorts out all the input spacepoints based on their distance to the z-axis
/// 2. Create seeds from 3 spacepoints with a small deviation from a striagh line
/// 3. Returns vertex position as (0,0,Z), where "Z" is the peak where seeds point to
template <typename spacepoint_t>
class ZScanSeedVertexFinder {
    public:
        /// Configuration struct
        struct Config{
            // maximum deviation in z-r plane between the first 2 spacepoints and the last 2 spacepoints
            float maxZRdeviation = 0.1;
            // maximum deviation in z-r plane between the first 2 spacepoints and the last 2 spacepoints
            float maxXYdeviation = 0.1;

            // thresholds for near, middle, and far spacepoints
            float rMinMiddle = 60.f * Acts::UnitConstants::mm;
            float rMaxMiddle = 120.f * Acts::UnitConstants::mm;

            // minimum slope in z-r plane for near-middle spacepoints, effectively removing triplets with large |eta|
            float minZRslope = 0.2;

            // bin size for the histogram
            float zBinSize = 1.f * Acts::UnitConstants::mm;

            // maximum z position of the vertex
            float maxZPosition = 300.f * Acts::UnitConstants::mm;
        };
        
        /// Const access to the config
        const Config& config() const { return m_cfg; }

        /// @brief Constructor
        ///
        /// @param cfg Configuration object 
        /// @param logger Logging instance
        ZScanSeedVertexFinder(const Config& cfg,
                              std::unique_ptr<const Logger> lgr = getDefaultLogger("ZScanSeedVertexFinder", Logging::DEBUG));

        /// @brief Destructor
        ~ZScanSeedVertexFinder() = default;

        /// @brief Finds the vertex based on the provided spacepoints
        /// @param spacepoints Vector of the input spacepoints; they do not need to be sorted anyhow
        /// @return Position of the vertex along z-axis
        std::vector<std::vector<float>> findVertex(const std::vector<spacepoint_t>& spacepoints) const;

    private:
        /// @brief Struct to store spacepoint combinations from near, middle, and far parts of the detector
        struct Triplet{
            Triplet(const spacepoint_t* aa, const spacepoint_t* bb, const spacepoint_t* cc) : a(aa), b(bb), c(cc) {}
            const spacepoint_t *a, *b, *c;
        };

        /// Configuration instance 
        Config m_cfg;

        /// @brief Sorts spacepoints into a separate vectors for near, middle, and far spacepoints
        /// @param spacepoints Vector of the input spacepoints;
        /// @return Vector of vectors for each set of spacepoints
        std::vector<std::vector<const spacepoint_t*>> sortSpacepoints(const std::vector<spacepoint_t>& spacepoints) const;

        /// @brief Makes triplets from the provided vectors of near, middle, and far spacepoints
        /// @param sorted_spacepoints Vector of input vector of vectors for each set of spacepoints
        /// @return Vector of valid triplets
        std::vector<Triplet> findTriplets(const std::vector<std::vector<const spacepoint_t*>>& sorted_spacepoints) const;

        /// @brief Validate the triplet based on "maxZRdeviation" and "maxXYdeviation"
        /// @param triplet A single triplet to be validated
        /// @return True if the deviations are within configured ranges
        bool isTripletValid(const Triplet triplet) const;

        /// @brief Calculates equation of the plane (alpha*x + beta*y + gamma*z + delta = 0), given the three points
        /// @param triplet A single triplet (with 3 spacepoints)
        /// @return array of {alpha,beta,gamma,delta}
        std::array<double,4> makePlaneFromTriplet(const Triplet triplet) const;

        /// @brief Find a point (=the vertex) that has minimum chi^2 with respect to all planes defined by the triplets
        /// @param triplets Vector of all valid triplets
        /// @return Position {x,y,z} of the vertex
        std::vector<float> findClosestPointFromPlanes(const std::vector<Triplet>& triplets) const;

        /// @brief Calculates parameters of the line (starting point + direction), given the three points
        /// @param triplet A single triplet (with 3 spacepoints)
        /// @return array of {starting_point_X,starting_point_Y,starting_point_Z, direction_X,direction_Y,direction_Z}
        std::array<double,6> makeLineFromTriplet(const Triplet triplet) const;

        /// @brief Find a point (=the vertex) that has minimum chi^2 with respect to all lines fitted through the triplets
        /// @param triplets Vector of all valid triplets
        /// @return Position {x,y,z} of the vertex
        std::vector<float> findClosestPointFromLines(const std::vector<Triplet>& triplets) const;

        std::array<double,2> makeZRLineFromTriplet(const Triplet triplet) const;
        std::vector<std::vector<float>> findClosestPointFromZRLines(const std::vector<Triplet>& triplets) const;

        /// @brief Creates a vector pretending to be a histogram of the estimated origins of the triplets
        /// @param triplets Vector of all valid triplets
        /// @return Histogram of the estimated origins
        std::vector<int> makeZHist(const std::vector<Triplet>& triplets) const;
        
        /// @brief Finds a single peak in the histogram
        /// @param hist Histograms with (preferenbly) a single peak
        /// @return Position along z-axis of the peak, taking into account "zBinSize"
        float findZPeak(const std::vector<int>& hist, int doSecondPeak) const;

        /// Logging instance
        std::unique_ptr<const Logger> m_logger;

        /// Private access to logging instance
        const Logger& logger() const { return *m_logger; }
};

} // namespace Acts

#include "Acts/Vertexing/ZScanSeedVertexFinder.ipp"
