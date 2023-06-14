// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Ray.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

namespace Acts {
/// @class SeedVertexFinder
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
class SeedVertexFinder {
 public:
  /// Configuration struct
  struct Config {
    // maximum deviation in phi between the near and middle spacepoints or
    // middle and far spacepoints
    float maxPhideviation = 0.08;
    // maximum deviation in X-Y between the first 2 spacepoints and the last 2
    // spacepoints
    float maxXYdeviation = 0.08;
    // maximum deviation in 3D between the first 2 spacepoints and the last 2
    // spacepoints
    float maxXYZdeviation = 0.08;

    // minimum angle between Z axis and a triplet, effectively removing triplets
    // with large |eta|
    float minTheta = 1.;

    // thresholds for near, middle, and far spacepoints
    float rMinNear = 20.f * Acts::UnitConstants::mm;
    float rMaxNear = 60.f * Acts::UnitConstants::mm;
    float rMinMiddle = 150.f * Acts::UnitConstants::mm;
    float rMaxMiddle = 190.f * Acts::UnitConstants::mm;
    float rMinFar = 280.f * Acts::UnitConstants::mm;
    float rMaxFar = 320.f * Acts::UnitConstants::mm;

    // number of phi slices, at least 3 to avoid duplicities
    // it should be less than 2*pi/maxPhideviation in order not to loop over
    // triplets that will be rejected by maxPhideviation anyway
    int numPhiSlices = 60;
    // use only a fraction of available phi slices to speed up calculations;
    float useFracPhiSlices = 0.5;

    // number of z slices
    int numZSlices = 150;
    // use only a fraction of available z slices to speed up calculations;
    float useFracZSlices = 0.5;
    // maximum |z| to consider, z slices will be done within the range
    // (-maxAbsZ,maxAbsZ) values of maxAbsZ, maxZPosition, rMaxFar, and minTheta
    // should be set reasonably with respect to each other
    float maxAbsZ = 450. * Acts::UnitConstants::mm;

    // maximum Z position of the vertex at the point closest to the Z axis
    float maxZPosition = 200.f * Acts::UnitConstants::mm;
    // maximum R position of the vertex at the point closest to the Z axis
    float maxRPosition = 10.f * Acts::UnitConstants::mm;

    // chi^2 minimalization will happen with respect to "planes" or "rays"
    std::string minimalizeWRT = "planes";

    // maximum number of iterations when discarding triplets with the largest
    // chi^2
    int maxIterations = 20;
    // each iteration, discard this fraction of triplets with the largest chi^2
    float removeFraction = 0.10;
    // if the vertex estimation moves less than this, stop iterations
    float minVtxShift = 0.3f * Acts::UnitConstants::mm;
  };

  /// Const access to the config
  const Config& config() const { return m_cfg; }

  /// @brief Constructor
  ///
  /// @param cfg Configuration object
  /// @param lgr Logging instance
  SeedVertexFinder(const Config& cfg,
                   std::unique_ptr<const Logger> lgr =
                       getDefaultLogger("SeedVertexFinder", Logging::DEBUG));

  /// @brief Destructor
  ~SeedVertexFinder() = default;

  /// @brief Finds the vertex based on the provided spacepoints
  /// @param spacepoints Vector of the input spacepoints; they do not need to be sorted anyhow
  /// @return Position of the vertex along z-axis
  Acts::Vector3 findVertex(const std::vector<spacepoint_t>& spacepoints);

 private:
  /// @brief Struct to store spacepoint combinations from near, middle, and far parts of the detector
  struct Triplet {
    Triplet(const spacepoint_t& aa, const spacepoint_t& bb,
            const spacepoint_t& cc)
        : a(aa),
          b(bb),
          c(cc),
          ray(Acts::Vector3::Zero(), Acts::Vector3::Zero()) {}
    const spacepoint_t &a, &b, &c;

    Acts::Ray3D ray;
  };

  /// Configuration instance
  Config m_cfg;

  /// @brief Sorts spacepoints into a separate vectors for near, middle, and far spacepoints; and for each slice of phi
  /// @param spacepoints Vector of the input spacepoints;
  /// @return Vector of vectors for each set of spacepoints
  std::vector<std::vector<std::vector<
      std::vector<std::pair<spacepoint_t const*, Acts::ActsScalar>>>>>
  sortSpacepoints(const std::vector<spacepoint_t>& spacepoints) const;

  /// @brief
  /// @param a
  /// @param b
  /// @return
  static bool compareSpacepoints(
      std::pair<spacepoint_t const*, Acts::ActsScalar>& a,
      std::pair<spacepoint_t const*, Acts::ActsScalar>& b);

  /// @brief Makes triplets from the provided vectors of near, middle, and far spacepoints; and for each slice of phi
  /// @param sorted_spacepoints Vector of input vector of vectors for each set of spacepoints
  /// @return Vector of valid triplets
  std::vector<Triplet> findTriplets(
      const std::vector<std::vector<std::vector<
          std::vector<std::pair<spacepoint_t const*, Acts::ActsScalar>>>>>&
          sorted_spacepoints);

  /// @brief Validate the triplet based on "maxXYdeviation", "maxXYZdeviation", "maxZPosition", and "maxRPosition"
  /// @param triplet A single triplet to be validated
  /// @return True if the deviations and fitted ray are within the configured ranges
  ///         If "minimalizeWRT"=="rays", then the fitted ray is also saved to
  ///         the triplet for later
  bool isTripletValid(Triplet& triplet);

  /// @brief Calculates equation of the plane (alpha*x + beta*y + gamma*z + delta = 0), given the three points
  /// @param triplet A single triplet (with 3 spacepoints)
  /// @return A pair of {{alpha,beta,gamma},delta}
  std::pair<Acts::Vector3, Acts::ActsScalar> makePlaneFromTriplet(
      const Triplet triplet) const;

  /// @brief Find a point (=the vertex) that has minimum chi^2 with respect to all planes defined by the triplets
  /// @param triplets Vector of all valid triplets
  /// @return Position {x,y,z} of the vertex
  Acts::Vector3 findClosestPointFromPlanes(
      const std::vector<Triplet>& triplets) const;

  /// @brief Calculates parameters of the ray (starting point + direction), given the three points
  /// @param triplet A single triplet (with 3 spacepoints)
  /// @return A ray of {starting_point, direction}
  Acts::Ray3D makeRayFromTriplet(const Triplet triplet) const;

  /// @brief Find a point (=the vertex) that has minimum chi^2 with respect to all rays fitted through the triplets
  /// @param triplets Vector of all valid triplets
  /// @return Position {x,y,z} of the vertex
  Acts::Vector3 findClosestPointFromRays(
      const std::vector<Triplet>& triplets) const;

  /// Logging instance
  std::unique_ptr<const Logger> m_logger;

  /// Private access to logging instance
  const Logger& logger() const { return *m_logger; }
};

}  // namespace Acts

#include "Acts/Vertexing/SeedVertexFinder.ipp"
