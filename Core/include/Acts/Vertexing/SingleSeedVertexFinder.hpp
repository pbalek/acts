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
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <algorithm>
#include <functional>
#include <string>
#include <utility>
#include <vector>

namespace Acts {

/// @class SingleSeedVertexFinder
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
class SingleSeedVertexFinder {
 public:
  /// Configuration struct
  struct Config {

    std::uint32_t minNumTriplets = 100;

    /// maximum deviation in phi between the near and middle spacepoints or
    /// middle and far spacepoints
    Acts::ActsScalar maxPhideviation = 0.08;
    /// maximum deviation in X-Y between the first 2 spacepoints and the last 2
    /// spacepoints
    Acts::ActsScalar maxXYdeviation = 0.08;
    /// maximum deviation in 3D between the first 2 spacepoints and the last 2
    /// spacepoints
    Acts::ActsScalar maxXYZdeviation = 0.08;

    /// minimum angle between Z axis and a triplet, effectively removing
    /// triplets with large |eta|
    Acts::ActsScalar minTheta = 1.;

    /// thresholds for near, middle, and far spacepoints
    Acts::ActsScalar rMinNear = 20.f * Acts::UnitConstants::mm;
    Acts::ActsScalar rMaxNear = 60.f * Acts::UnitConstants::mm;
    Acts::ActsScalar rMinMiddle = 150.f * Acts::UnitConstants::mm;
    Acts::ActsScalar rMaxMiddle = 190.f * Acts::UnitConstants::mm;
    Acts::ActsScalar rMinFar = 280.f * Acts::UnitConstants::mm;
    Acts::ActsScalar rMaxFar = 320.f * Acts::UnitConstants::mm;

    /// number of phi slices, at least 3 to avoid duplicities
    /// it should be less than 2*pi/maxPhideviation in order not to loop over
    /// triplets that will be rejected by maxPhideviation anyway
    std::uint32_t numPhiSlices = 60;
    /// use only a fraction of available phi slices to speed up calculations;
    Acts::ActsScalar useFracPhiSlices = 1.0;

    /// number of z slices
    std::uint32_t numZSlices = 150;
    /// use only a fraction of available z slices to speed up calculations;
    Acts::ActsScalar useFracZSlices = 1.0;
    /// maximum |z| to consider, z slices will be done within the range
    /// (-maxAbsZ,maxAbsZ)
    /// values of maxAbsZ, maxZPosition, rMaxFar, and minTheta should be
    /// set reasonably with respect to each other
    Acts::ActsScalar maxAbsZ = 450. * Acts::UnitConstants::mm;

    /// maximum Z position of the vertex at the point closest to the Z axis
    Acts::ActsScalar maxZPosition = 200.f * Acts::UnitConstants::mm;
    /// maximum R position of the vertex at the point closest to the Z axis
    Acts::ActsScalar maxRPosition = 10.f * Acts::UnitConstants::mm;

    /// chi^2 minimalization will happen with respect to "planes" or "rays" or "mixed"
    std::string minimalizeWRT = "mixed";

    /// maximum number of iterations when discarding triplets with the largest
    /// chi^2
    std::uint32_t maxIterations = 30;
    /// minimum number of iterations when discarding triplets with the largest
    /// chi^2
    std::uint32_t minIterations = 10;
    /// each iteration, discard this fraction of triplets with the largest chi^2
    Acts::ActsScalar removeFraction = 0.1;
    /// if the vertex estimation moves less than this, stop iterations
    Acts::ActsScalar minVtxShift = 0.3f * Acts::UnitConstants::mm;
    /// valid only for "mixed" option of "minimalizeWRT"
    /// eccentricity of a virtual ellipse: semi-minor axis defines the
    /// unity distance of chi2 off-plane and semi-major axis defines the
    /// unity distance of chi2 in-plane
    Acts::ActsScalar mixedEccentricity = 0.75;
  };

  /// Const access to the config
  const Config& config() const { return m_cfg; }

  /// @brief Constructor
  /// @param cfg Configuration object
  /// @param lgr Logging instance
  SingleSeedVertexFinder(const Config& cfg,
                         std::unique_ptr<const Logger> lgr = getDefaultLogger(
                             "SingleSeedVertexFinder", Logging::INFO));

  /// @brief Destructor
  ~SingleSeedVertexFinder() = default;

  /// @brief Finds the vertex based on the provided spacepoints
  /// @param spacepoints Vector of the input spacepoints; they do not need to be sorted anyhow
  /// @return Position of the vertex
  // Acts::Result<Acts::Vector3> findVertex(
  std::pair<Acts::Result<Acts::Vector3>, std::vector<std::vector<Acts::ActsScalar>>> findVertex(
      const std::vector<spacepoint_t>& spacepoints) const;

 private:
  /// @brief Struct to store spacepoint combinations from near, middle, and far parts of the detector. Also stores straight line fit through the spacepoints in case minimalizeWRT=="rays", so it's not fitted twice
  struct Triplet {
    Triplet(Acts::Vector3 abg, Acts::ActsScalar delta, Acts::Vector3 start, Acts::Vector3 dir, 
            const spacepoint_t& aa, const spacepoint_t& bb, const spacepoint_t& cc)
        : planeABG(abg), planeDelta(delta), startPoint(start), direction(dir), 
          distance(0.), nearSP(aa), middleSP(bb), farSP(cc) {}
    Triplet(const spacepoint_t& aa, const spacepoint_t& bb, const spacepoint_t& cc)
        : distance(-2.), nearSP(aa), middleSP(bb), farSP(cc) {} // TODO: unnecessary assignment
    Triplet(Acts::ActsScalar dist, const spacepoint_t& aa, const spacepoint_t& bb, const spacepoint_t& cc)
        : distance(dist), nearSP(aa), middleSP(bb), farSP(cc) {}

    Acts::Vector3 getStartPoint() const {
      return startPoint;
    }
    Acts::Vector3 getDirection() const {
      return direction;
    }

    Acts::Vector3 getPlaneABG() const {
      return planeABG;
    }
    Acts::ActsScalar getPlaneDelta() const {
      return planeDelta;
    }

    Acts::ActsScalar getDistance() const {
      return distance;
    }
    void setDistance(Acts::ActsScalar d) {
      distance=d;
    }
    
    Acts::Vector3 planeABG;
    Acts::ActsScalar planeDelta;
    Acts::Vector3 startPoint, direction;
    Acts::ActsScalar distance;
    std::reference_wrapper<const spacepoint_t> nearSP, middleSP, farSP;
  };

  /// @brief Struct to store sorted spacepoints for each layer (near, middle, and far), for each slice of phi, and for each slice of z
  struct SortedSpacepoints {
    SortedSpacepoints(const int phiSize, const int zSize) {
      std::vector<std::pair<spacepoint_t const*, Acts::ActsScalar>> helper = {};
      std::vector<std::vector<std::pair<spacepoint_t const*, Acts::ActsScalar>>>
          helperZ(zSize, helper);
      std::vector<std::vector<
          std::vector<std::pair<spacepoint_t const*, Acts::ActsScalar>>>>
          helperPhi(phiSize, helperZ);
      sortedSP.fill(helperPhi);
    }

    /// @brief Provides non-const vector of spacepoints for a given layer, phi slice, and z slice
    /// @param layer Index of the layer (near=0, middle=1, far=2)
    /// @param phi Index of the phi slice
    /// @param z Index of the z slice
    /// @return Non-const vector of spacepoints
    inline std::vector<std::pair<spacepoint_t const*, Acts::ActsScalar>>& addSP(
        int layer, int phi, int z) {
      return sortedSP[layer][phi][z];
    }

    /// @brief Provides const vector of spacepoints for a given layer, phi slice, and z slice
    /// @param layer Index of the layer (near=0, middle=1, far=2)
    /// @param phi Index of the phi slice
    /// @param z Index of the z slice
    /// @return Const vector of spacepoints
    inline const std::vector<std::pair<spacepoint_t const*, Acts::ActsScalar>>&
    getSP(int layer, int phi, int z) const {
      return sortedSP[layer][phi][z];
    }

    std::array<std::vector<std::vector<std::vector<
                   std::pair<spacepoint_t const*, Acts::ActsScalar>>>>,
               3>
        sortedSP;
  };

  /// Configuration instance
  Config m_cfg;

  /// @brief Sorts spacepoints into a separate vectors for near, middle, and far spacepoints; for each slice of phi; and for each slice of z
  /// @param spacepoints Vector of the input spacepoints;
  /// @return Struct of the sorted spacepoints
  Acts::SingleSeedVertexFinder<spacepoint_t>::SortedSpacepoints sortSpacepoints(
      const std::vector<spacepoint_t>& spacepoints) const;

  /// @brief Makes triplets from the provided vectors of near, middle, and far spacepoints; for each slice of phi; and for each slice of z
  /// @param sortedSpacepoints Struct of the sorted spacepointss
  /// @return Vector of valid triplets
  std::vector<Triplet> findTriplets(
      const Acts::SingleSeedVertexFinder<spacepoint_t>::SortedSpacepoints&
          sortedSpacepoints) const;

  /// @brief Validate the triplet based on "maxXYdeviation", "maxXYZdeviation", "maxZPosition", and "maxRPosition"
  /// @param a,b,c Spacepoints to form a triplet to be validated
  /// @return True if the deviations and fitted ray are within the configured ranges
  ///         If "minimalizeWRT"!="planes", then the fitted ray is also saved to
  ///         the triplet for later; if "minimalizeWRT"!="rays", then the fitted plane
  ///         is also saved for later
  Triplet tripletValidationAndFit(const spacepoint_t& a, const spacepoint_t& b, const spacepoint_t& c) const;

  /// @brief Fit the triplet with a ray and a plane
  /// @param a,b,c Spacepoints to form a triplet
  /// @return A fitted triplet with "ray" and "plane" saved for later;
  ///         One of these will be skipped if it's sure it will not be needed
  Triplet fitTriplet(
    const spacepoint_t& spA, const spacepoint_t& spB, const spacepoint_t& spC) const;

  /// @brief Calculates equation of the plane (alpha*x + beta*y + gamma*z + delta = 0), given the three points
  /// @param a,b,c Spacepoints to form a triplet to be validated
  /// @return A pair of {{alpha,beta,gamma},delta}
  std::pair<Acts::Vector3, Acts::ActsScalar> makePlaneFromTriplet(
      const spacepoint_t& a, const spacepoint_t& b, const spacepoint_t& c) const;

  // /// @brief Find a point (=the vertex) that has minimum chi^2 with respect to all planes defined by the triplets
  // /// @param triplets Vector of all valid triplets
  // /// @return Position {x,y,z} of the vertex
  // Acts::Vector3 findClosestPointFromPlanes(
  //     std::vector<typename Acts::SingleSeedVertexFinder<spacepoint_t>::Triplet>& triplets) const;

  /// @brief Calculates parameters of the ray (starting point + direction), given the three points
  /// @param a,b,c Spacepoints to form a triplet to be validated
  /// @return A ray of {starting_point, direction}
  std::pair<Acts::Vector3, Acts::Vector3> makeRayFromTriplet(const spacepoint_t& a, const spacepoint_t& b, const spacepoint_t& c) const;

  // /// @brief Find a point (=the vertex) that has minimum chi^2 with respect to all rays fitted through the triplets
  // /// @param triplets Vector of all valid triplets
  // /// @return Position {x,y,z} of the vertex
  // Acts::Vector3 findClosestPointFromRays(
  //     std::vector<typename Acts::SingleSeedVertexFinder<spacepoint_t>::Triplet>& triplets) const;

  /// @brief Find a point (=the vertex) that has minimum chi^2 with respect to all rays/planes fitted through the triplets, while accounting for mixedEccentricity
  /// @param triplets Vector of all valid triplets
  /// @return Position {x,y,z} of the vertex
  Acts::Vector3 findClosestPoint(
    std::vector<typename Acts::SingleSeedVertexFinder<spacepoint_t>::Triplet>& allTriples, std::vector<std::vector<Acts::ActsScalar>>& rejectVector, std::vector<std::vector<Acts::ActsScalar>>& vtx_iter) const;

  /// @brief Square of effective eccentricity; it is set to the value of m_cfg.mixedEccentricity^2 for minimalizeWRT="mixed", and to the values of 0 and 1 for "rays" and "planes", respectively.
  Acts::ActsScalar m_effectEccSq;

  /// Logging instance
  std::unique_ptr<const Logger> m_logger;

  /// Private access to logging instance
  const Logger& logger() const { return *m_logger; }
};

}  // namespace Acts

#include "Acts/Vertexing/SingleSeedVertexFinder.ipp"
