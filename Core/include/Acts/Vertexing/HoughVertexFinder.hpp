// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding/HoughTransformUtils.hpp"
#include "Acts/Seeding/HoughVectors.hpp"
#include "Acts/Utilities/AlgebraHelpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Ray.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <algorithm>
#include <numeric>
#include <functional>
#include <string>
#include <utility>
#include <vector>
#include <execution>

namespace Acts {

/// @class HoughVertexFinder
///
/// @brief Implements the vertex finder based on the track seeds
/// 0. Assumes there is only 1 vertex and that it has a high multiplicity
/// 1. Sorts out all the input spacepoints based on their distance to the z-axis
/// 2. Create seeds from 3 spacepoints with a small deviation from a straigh
/// line
/// 3. Find a point with a minimal distance from either planes
/// (minimalizeWRT="planes") or rays (minimalizeWRT="rays") defined by the the
/// seeds
/// 4. Returns the point position as the vertex
template <typename spacepoint_t>
class HoughVertexFinder {
 public:
  /// Configuration struct
  struct Config {

    /// Ideal amount of spacepoints; |eta| range will be limited to
    /// contain approximately this amount of SPs
    std::uint32_t targetSPs = 10000;

    /// Minimum and maximum ranges in |eta|; the |eta| will no be 
    /// set outside these bounds even if targetSPs is not reached
    Acts::ActsScalar minAbsEta = 0.3f;
    Acts::ActsScalar maxAbsEta = 3.0f;

    /// Minimum number of hits in Hough plane to consider
    /// the cell to contain a track
    std::uint32_t minHits = 4;

    /// Number of neighbouring bins in Hough plane
    /// to fill in the cot(theta) direction
    std::uint32_t fillNeighbours = 0;

    /// Approximate amount of measurements in the |eta| range expressed 
    /// as a fraction of all measurements within the whole |eta| range.
    /// Measurements are assumed to be distributed uniformly with 
    /// the |eta| range. The first |eta| range starts at 0., the others 
    /// start at the endpoint of the previous range.
    std::vector<Acts::ActsScalar> absEtaRanges{2.,4.};
    std::vector<Acts::ActsScalar> absEtaFractions{0.4,0.6};

    /// Iterations along Z axis. The algorithm may change both the range 
    /// in Z (to reduce time) and the number of bins in Z (to achieve
    /// more precise result).
    std::vector<Acts::ActsScalar> rangeIterZ{200.f * Acts::UnitConstants::mm, 
                                              30.f * Acts::UnitConstants::mm, 
                                              16.f * Acts::UnitConstants::mm};
    std::vector<std::uint32_t> nBinsZIterZ{800, 180, 80};
    std::vector<std::uint32_t> nBinsCotThetaIterZ{8000, 8000, 8000};

    /// For every magnitude (in natural logarithm) below targetSPs, the number of
    /// bins in cot(theta) will decrease by this factor. Thus, the actual number
    /// of bins can be smaller than stated in "nBinsCotThetaIterZ".
    Acts::ActsScalar binsCotThetaDecrease{1.35};

    /// Width of the peak used for more precise z-position estimate
    std::uint32_t peakWidth{3};

    /// Default position of the vertex in X, Y, and Z coordinates
    Acts::Vector3 defVtxPosition{0.f * Acts::UnitConstants::mm,
                                 0.f * Acts::UnitConstants::mm,
                                 0.f * Acts::UnitConstants::mm};
  };

  /// Const access to the config
  const Config& config() const { return m_cfg; }

  /// @brief Constructor
  /// @param cfg Configuration object
  /// @param lgr Logging instance
  HoughVertexFinder(const Config& cfg,
                         std::unique_ptr<const Logger> lgr = getDefaultLogger(
                             "HoughVertexFinder", Logging::INFO));

  /// @brief Destructor
  ~HoughVertexFinder() = default;

  /// @brief Finds the vertex based on the provided spacepoints
  /// @param spacepoints Vector of the input spacepoints; they do not need to be sorted anyhow
  /// @return Position of the vertex
  Acts::Result<Acts::Vector3> find(
      const std::vector<spacepoint_t>& spacepoints) const;

 private:

  /// Configuration instance
  Config m_cfg;

  /// @brief Returns the positions of the peak along Z axis in the pojection of the Hough plane
  /// @param spacepoints Set of all spacepoints within the event
  /// @param vtxOld Previous position of the vertex
  /// @param rangeZ Range in along Z around vtxOld_z to consider when looking for the new vertex
  /// @param numZBins Number of bins along Z axis
  /// @param minCotTheta Minimum theta to consider for the spacepoint
  /// @param maxCotTheta Maximum theta to consider for the spacepoint
  /// @param numCotThetaBins Number of bins along cot(theta) axis
  /// @return Position of the vertex in (X,Y,Z)
  Acts::Result<Acts::Vector3> findHoughPeak(const std::vector<spacepoint_t>& spacepoints,
    Acts::Vector3 vtxOld, Acts::ActsScalar rangeZ, std::uint32_t numZBins,
    Acts::ActsScalar minCotTheta, Acts::ActsScalar maxCotTheta, std::uint32_t numCotThetaBins) const;

  /// Logging instance
  std::unique_ptr<const Logger> m_logger;

  /// Private access to logging instance
  const Logger& logger() const { return *m_logger; }
};

}  // namespace Acts

#include "Acts/Vertexing/HoughVertexFinder.ipp"
