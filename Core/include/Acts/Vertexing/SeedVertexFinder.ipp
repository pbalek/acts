// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cmath>

#include <Eigen/Eigenvalues>

template <typename spacepoint_t>
Acts::SeedVertexFinder<spacepoint_t>::SeedVertexFinder(
    const Acts::SeedVertexFinder<spacepoint_t>::Config& cfg,
    std::unique_ptr<const Logger> lgr)
    : m_cfg(cfg), m_logger(std::move(lgr)) {
  if (std::isnan(cfg.maxPhideviation))
    ACTS_ERROR("value of maxPhideviation was not initialized");
  if (std::isnan(cfg.maxXYdeviation))
    ACTS_ERROR("value of maxXYdeviation was not initialized");
  if (std::isnan(cfg.maxXYZdeviation))
    ACTS_ERROR("value of maxXYZdeviation was not initialized");
  if (std::isnan(cfg.minTheta))
    ACTS_ERROR("value of minTheta was not initialized");
  if (std::isnan(cfg.rMinNear))
    ACTS_ERROR("value of rMinNear was not initialized");
  if (std::isnan(cfg.rMaxNear))
    ACTS_ERROR("value of rMaxNear was not initialized");
  if (std::isnan(cfg.rMinMiddle))
    ACTS_ERROR("value of rMinMiddle was not initialized");
  if (std::isnan(cfg.rMaxMiddle))
    ACTS_ERROR("value of rMaxMiddle was not initialized");
  if (std::isnan(cfg.rMinFar))
    ACTS_ERROR("value of rMinFar was not initialized");
  if (std::isnan(cfg.rMaxFar))
    ACTS_ERROR("value of rMaxFar was not initialized");
  if (std::isnan(cfg.numPhiSlices))
    ACTS_ERROR("value of numPhiSlices was not initialized");
  if (cfg.numPhiSlices < 3)
    ACTS_INFO("value of numPhiSlices is "
              << cfg.numPhiSlices
              << ", which is less than 3. There will be duplicate triplets.");
  if (std::isnan(cfg.useFracPhiSlices))
    ACTS_ERROR("value of useFracPhiSlices was not initialized");
  if (cfg.useFracPhiSlices < 0. || cfg.useFracPhiSlices >= 1.)
    ACTS_ERROR("value of useFracPhiSlices is "
               << cfg.useFracPhiSlices
               << ", allowed values are between 0 and 1");
  if (std::isnan(cfg.numZSlices))
    ACTS_ERROR("value of numZSlices was not initialized");
  if (std::isnan(cfg.useFracZSlices))
    ACTS_ERROR("value of useFracZSlices was not initialized");
  if (cfg.useFracZSlices < 0. || cfg.useFracZSlices >= 1.)
    ACTS_ERROR("value of useFracZSlices is "
               << cfg.useFracZSlices << ", allowed values are between 0 and 1");
  if (std::isnan(cfg.maxAbsZ))
    ACTS_ERROR("value of maxAbsZ was not initialized");
  if (std::isnan(cfg.maxZPosition))
    ACTS_ERROR("value of maxZPosition was not initialized");
  if (std::isnan(cfg.maxRPosition))
    ACTS_ERROR("value of maxRPosition was not initialized");
  if (cfg.minimalizeWRT != "planes" && cfg.minimalizeWRT != "rays") {
    ACTS_ERROR("value of minimalizeWRT is "
               << cfg.minimalizeWRT
               << ", allowed values are \"planes\" or \"rays\" ");
  }
  if (std::isnan(cfg.maxIterations))
    ACTS_ERROR("value of maxIterations was not initialized");
  if (std::isnan(cfg.removeFraction))
    ACTS_ERROR("value of removeFraction was not initialized");
  if (cfg.removeFraction < 0. || cfg.removeFraction >= 1.)
    ACTS_ERROR("value of removeFraction is "
               << cfg.removeFraction << ", allowed values are between 0 and 1");
  if (std::isnan(cfg.minVtxShift))
    ACTS_ERROR("value of minVtxShift was not initialized");
}

template <typename spacepoint_t>
Acts::Vector3 Acts::SeedVertexFinder<spacepoint_t>::findVertex(
    const std::vector<spacepoint_t>& spacepoints) {
  // sort spacepoints to different phi and z bins
  std::vector<std::vector<std::vector<
      std::vector<std::pair<spacepoint_t const*, Acts::ActsScalar>>>>>
      sorted_spacepoints = sortSpacepoints(spacepoints);

  // find triplets
  std::vector<Acts::SeedVertexFinder<spacepoint_t>::Triplet> triplets =
      findTriplets(sorted_spacepoints);

  // if no valid triplets found
  if (triplets.empty())
    return {};

  Acts::Vector3 vtx = Acts::Vector3::Zero();
  if (m_cfg.minimalizeWRT == "planes") {
    // find a point closest to all planes defined by the triplets
    vtx = findClosestPointFromPlanes(triplets);
  } else if (m_cfg.minimalizeWRT == "rays") {
    // find a point closest to all rays fitted through the triplets
    vtx = findClosestPointFromRays(triplets);
  } else {
    ACTS_ERROR("value of minimalizeWRT is "
               << m_cfg.minimalizeWRT
               << ", allowed values are \"planes\" or \"rays\" ");
  }

  return vtx;
}

template <typename spacepoint_t>
std::vector<std::vector<
    std::vector<std::vector<std::pair<spacepoint_t const*, Acts::ActsScalar>>>>>
Acts::SeedVertexFinder<spacepoint_t>::sortSpacepoints(
    const std::vector<spacepoint_t>& spacepoints) const {
  std::vector<std::pair<spacepoint_t const*, Acts::ActsScalar>> helper = {};
  std::vector<std::vector<std::pair<spacepoint_t const*, Acts::ActsScalar>>>
      helper_eta(m_cfg.numZSlices, helper);

  std::vector<std::vector<
      std::vector<std::pair<spacepoint_t const*, Acts::ActsScalar>>>>
      near_spacepoints(m_cfg.numPhiSlices, helper_eta),
      middle_spacepoints(m_cfg.numPhiSlices, helper_eta),
      far_spacepoints(m_cfg.numPhiSlices, helper_eta);

  for (const auto& sp : spacepoints) {
    // phi will be saved for later
    Acts::ActsScalar phi = detail::radian_pos(std::atan2(sp.y(), sp.x()));
    int phislice = (int)(phi / (2 * M_PI) * m_cfg.numPhiSlices);

    if (fabs(sp.z()) >= m_cfg.maxAbsZ)
      continue;
    int zslice = (int)((sp.z() + m_cfg.maxAbsZ) / (2 * m_cfg.maxAbsZ) *
                       m_cfg.numZSlices);

    // input spacepoint is sorted into one subset
    if (sp.r() < m_cfg.rMinMiddle) {
      if (m_cfg.rMinNear < sp.r() && sp.r() < m_cfg.rMaxNear) {
        if (std::fmod(m_cfg.useFracPhiSlices * phislice, 1.0) >=
            m_cfg.useFracPhiSlices)
          continue;
        near_spacepoints.at(phislice).at(zslice).emplace_back(
            (spacepoint_t const*)&sp, phi);
      }
    } else {
      if (sp.r() < m_cfg.rMinFar) {
        if (sp.r() < m_cfg.rMaxMiddle) {
          if (std::fmod(m_cfg.useFracZSlices * zslice, 1.0) >=
              m_cfg.useFracZSlices)
            continue;
          middle_spacepoints.at(phislice).at(zslice).emplace_back(
              (spacepoint_t const*)&sp, phi);
        }
      } else {
        if (sp.r() < m_cfg.rMaxFar) {
          far_spacepoints.at(phislice).at(zslice).emplace_back(
              (spacepoint_t const*)&sp, phi);
        }
      }
    }
  }

  return {near_spacepoints, middle_spacepoints, far_spacepoints};
}

template <typename spacepoint_t>
std::vector<typename Acts::SeedVertexFinder<spacepoint_t>::Triplet>
Acts::SeedVertexFinder<spacepoint_t>::findTriplets(
    const std::vector<std::vector<std::vector<
        std::vector<std::pair<spacepoint_t const*, Acts::ActsScalar>>>>>&
        sorted_spacepoints) {
  std::vector<Acts::SeedVertexFinder<spacepoint_t>::Triplet> triplets;

  int phistep =
      (int)(m_cfg.maxPhideviation / (2 * M_PI / m_cfg.numPhiSlices)) + 1;

  // calculate limits for middle spacepoints
  Acts::Vector2 vecA{-m_cfg.maxAbsZ + m_cfg.maxZPosition, m_cfg.rMinFar};
  vecA /= 2.;
  Acts::Vector2 vecB = {vecA[1], -vecA[0]};
  vecB /= std::tan(m_cfg.maxXYZdeviation);
  Acts::Vector2 posR = Acts::Vector2(-m_cfg.maxZPosition, 0.) + vecA + vecB;
  Acts::ActsScalar R =
      std::sqrt(vecA.dot(vecA)) / std::sin(m_cfg.maxXYZdeviation);
  Acts::ActsScalar constB = -2. * posR[0];
  Acts::ActsScalar constC =
      posR[0] * posR[0] +
      (posR[1] - m_cfg.rMaxNear) * (posR[1] - m_cfg.rMaxNear) - R * R;
  Acts::ActsScalar maxZMiddle =
      -1. * (-constB - sqrt(constB * constB - 4. * constC)) / 2.;

  // save some constant values for later
  Acts::ActsScalar rNearRatio[2] = {m_cfg.rMinNear / m_cfg.rMaxMiddle,
                                    m_cfg.rMaxNear / m_cfg.rMinMiddle};
  Acts::ActsScalar rMiddle[2] = {m_cfg.rMaxMiddle, m_cfg.rMinMiddle};
  Acts::ActsScalar rFarDelta[2] = {
      m_cfg.rMaxFar - m_cfg.rMinMiddle,
      m_cfg.rMinFar - m_cfg.rMaxMiddle,
  };
  Acts::ActsScalar zBinLength = 2. * m_cfg.maxAbsZ / m_cfg.numZSlices;

  // limits in terms of slice numbers
  int limitMiddleSlice_from = (int)((-maxZMiddle + m_cfg.maxAbsZ) / zBinLength);
  int limitMiddleSlice_to =
      (int)((maxZMiddle + m_cfg.maxAbsZ) / zBinLength + 1);
  int limitAbsZSlice_from =
      (int)((-m_cfg.maxZPosition + m_cfg.maxAbsZ) / zBinLength + 0.01);
  int limitAbsZSlice_to =
      (int)((m_cfg.maxZPosition + m_cfg.maxAbsZ) / zBinLength + 1.01);

  for (int middlez = limitMiddleSlice_from; middlez < limitMiddleSlice_to;
       ++middlez) {
    // skip slices that are empty anyway
    if (std::fmod(m_cfg.useFracZSlices * middlez, 1.0) >= m_cfg.useFracZSlices)
      continue;

    // calculate limits for near spacepoints, assuming the middle spacepoints
    // are within some boundaries
    bool isLess_from = (middlez < limitAbsZSlice_from);
    int nearz_from;
    float deltaZfrom = (middlez - limitAbsZSlice_from - 1) * zBinLength;
    float angleZfrom =
        std::atan2(rMiddle[isLess_from], deltaZfrom) + m_cfg.maxXYZdeviation;
    if (angleZfrom > M_PI)
      nearz_from = 0;
    else {
      float new_deltaZfrom =
          rMiddle[isLess_from] / std::tan(angleZfrom) / zBinLength;
      nearz_from =
          (int)((new_deltaZfrom)*rNearRatio[isLess_from]) + limitAbsZSlice_from;
      if (nearz_from < 0)
        nearz_from = 0;
    }

    bool isLess_to = (middlez < limitAbsZSlice_to);
    int nearz_to;
    float deltaZto = (middlez - limitAbsZSlice_to + 1) * zBinLength;
    float angleZto =
        std::atan2(rMiddle[!isLess_to], deltaZto) - m_cfg.maxXYZdeviation;
    if (angleZto < 0)
      nearz_to = m_cfg.numZSlices;
    else {
      float new_deltaZto =
          rMiddle[!isLess_to] / std::tan(angleZto) / zBinLength;
      nearz_to =
          (int)((new_deltaZto)*rNearRatio[!isLess_to]) + limitAbsZSlice_to;
      if (nearz_to > m_cfg.numZSlices)
        nearz_to = m_cfg.numZSlices;
    }

    for (int nearz = nearz_from; nearz < nearz_to; ++nearz) {
      // calculate limits for far spacepoits, assuming middle and near
      // spacepoits are within some boundaries
      bool isMiddleLess = (middlez < nearz);

      int farz_from;
      float delta2Zfrom = (middlez - nearz - 1) * zBinLength;
      float angle2Zfrom = std::atan2(rFarDelta[isMiddleLess], delta2Zfrom) +
                          m_cfg.maxXYZdeviation;
      if (angle2Zfrom > M_PI)
        farz_from = 0;
      else {
        farz_from = (int)(rFarDelta[isMiddleLess] / std::tan(angle2Zfrom) /
                          zBinLength) +
                    middlez;
        if (farz_from < 0)
          farz_from = 0;
        else if (farz_from >= m_cfg.numZSlices)
          continue;
      }

      int farz_to;
      float delta2Zto = (middlez - nearz + 1) * zBinLength;
      float angle2Zto = std::atan2(rFarDelta[!isMiddleLess], delta2Zto) -
                        m_cfg.maxXYZdeviation;
      if (angle2Zto < 0)
        farz_to = m_cfg.numZSlices;
      else {
        farz_to =
            (int)(rFarDelta[!isMiddleLess] / std::tan(angle2Zto) / zBinLength) +
            middlez + 1;
        if (farz_to > m_cfg.numZSlices)
          farz_to = m_cfg.numZSlices;
        else if (farz_to < 0)
          continue;
      }

      for (int farz = farz_from; farz < farz_to; farz++) {
        // loop over near phi slices
        for (int nearphi = 0; nearphi < m_cfg.numPhiSlices; ++nearphi) {
          // skip slices that are empty anyway
          if (std::fmod(m_cfg.useFracPhiSlices * nearphi, 1.0) >=
              m_cfg.useFracPhiSlices)
            continue;

          // loop over some middle phi slices
          for (int middlephi_h = nearphi - phistep;
               middlephi_h <= nearphi + phistep; ++middlephi_h) {
            int middlephi =
                (middlephi_h + m_cfg.numPhiSlices) % m_cfg.numPhiSlices;
            // loop over some far phi slices
            for (int farphi_h = middlephi - phistep;
                 farphi_h <= middlephi + phistep; ++farphi_h) {
              int farphi = (farphi_h + m_cfg.numPhiSlices) % m_cfg.numPhiSlices;

              // for all near spacepoints in this slice
              for (const auto& near_sp :
                   sorted_spacepoints.at(0).at(nearphi).at(nearz)) {
                Acts::ActsScalar phiA = near_sp.second;
                // for all middle spacepoints in this slice
                for (const auto& middle_sp :
                     sorted_spacepoints.at(1).at(middlephi).at(middlez)) {
                  Acts::ActsScalar phiB = middle_sp.second;
                  Acts::ActsScalar delta_phiAB =
                      detail::difference_periodic(phiA, phiB, 2 * M_PI);
                  if (std::abs(delta_phiAB) > m_cfg.maxPhideviation)
                    continue;
                  // for all far spacepoints in this slice
                  for (const auto& far_sp :
                       sorted_spacepoints.at(2).at(farphi).at(farz)) {
                    Acts::ActsScalar phiC = far_sp.second;
                    Acts::ActsScalar delta_phiBC =
                        detail::difference_periodic(phiB, phiC, 2 * M_PI);
                    if (std::abs(delta_phiBC) > m_cfg.maxPhideviation)
                      continue;

                    Acts::SeedVertexFinder<spacepoint_t>::Triplet tr(
                        *near_sp.first, *middle_sp.first, *far_sp.first);

                    if (isTripletValid(tr)) {
                      triplets.push_back(tr);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  return triplets;
}

template <typename spacepoint_t>
bool Acts::SeedVertexFinder<spacepoint_t>::isTripletValid(
    Acts::SeedVertexFinder<spacepoint_t>::Triplet& triplet) {
  // slope for near+middle spacepoints
  Acts::ActsScalar alpha1 =
      std::atan2(triplet.a.y() - triplet.b.y(), triplet.a.x() - triplet.b.x());
  // slope for middle+far spacepoints
  Acts::ActsScalar alpha2 =
      std::atan2(triplet.b.y() - triplet.c.y(), triplet.b.x() - triplet.c.x());
  // these two slopes shouldn't be too different
  Acts::ActsScalar delta_alpha =
      detail::difference_periodic(alpha1, alpha2, 2 * M_PI);
  if (std::abs(delta_alpha) > m_cfg.maxXYdeviation)
    return false;

  // near-middle ray
  Acts::Vector3 ab{triplet.a.x() - triplet.b.x(), triplet.a.y() - triplet.b.y(),
                   triplet.a.z() - triplet.b.z()};
  // middle-far ray
  Acts::Vector3 bc{triplet.b.x() - triplet.c.x(), triplet.b.y() - triplet.c.y(),
                   triplet.b.z() - triplet.c.z()};
  // dot product of these two
  Acts::ActsScalar costheta =
      (ab.dot(bc)) / (std::sqrt(ab.dot(ab)) * std::sqrt(bc.dot(bc)));
  Acts::ActsScalar theta = std::acos(costheta);
  if (theta > m_cfg.maxXYZdeviation)
    return false;

  // reject the ray if it doesn't come close to the z-axis
  Acts::Ray3D ray = makeRayFromTriplet(triplet);
  const Acts::Vector3& start_point = ray.origin();
  const Acts::Vector3& direction = ray.dir();

  Acts::ActsScalar tanTheta =
      std::sqrt(direction[0] * direction[0] + direction[1] * direction[1]) /
      direction[2];
  if (std::fabs(tanTheta) < std::tan(m_cfg.minTheta))
    return false;

  // norm to z-axis and to the ray
  Acts::Vector3 norm{-1. * direction[1], 1. * direction[0], 0};
  // nearest distance from the ray to z-axis
  Acts::ActsScalar dist =
      std::fabs(start_point.dot(norm)) / std::sqrt(norm.dot(norm));
  if (dist > m_cfg.maxRPosition)
    return false;

  // cross product of direction and norm
  Acts::Vector3 direction_x_n = direction.cross(norm);
  // z coordinate of the nearest distance from the ray to z-axis
  Acts::ActsScalar zdist = direction_x_n.dot(start_point) / (norm.dot(norm));
  if (std::fabs(zdist) > m_cfg.maxZPosition)
    return false;

  if (m_cfg.minimalizeWRT == "rays") {
    // save for later
    triplet.ray = ray;
  }

  return true;
}

template <typename spacepoint_t>
std::pair<Acts::Vector3, Acts::ActsScalar>
Acts::SeedVertexFinder<spacepoint_t>::makePlaneFromTriplet(
    const Acts::SeedVertexFinder<spacepoint_t>::Triplet triplet) const {
  Acts::Vector3 a{triplet.a.x(), triplet.a.y(), triplet.a.z()};
  Acts::Vector3 b{triplet.b.x(), triplet.b.y(), triplet.b.z()};
  Acts::Vector3 c{triplet.c.x(), triplet.c.y(), triplet.c.z()};

  Acts::Vector3 ba = b - a, ca = c - a;

  Acts::Vector3 abg = ba.cross(ca);
  // vector (alpha,beta,gamma) normalized to unity for convenience
  abg /= std::sqrt(abg.dot(abg));
  Acts::ActsScalar delta = -1. * abg.dot(a);

  // plane (alpha*x + beta*y + gamma*z + delta = 0), splitted to {{alpha, beta,
  // gamma}, delta} for convenience
  return {abg, delta};
}

template <typename spacepoint_t>
Acts::Vector3 Acts::SeedVertexFinder<spacepoint_t>::findClosestPointFromPlanes(
    const std::vector<Acts::SeedVertexFinder<spacepoint_t>::Triplet>& triplets)
    const {
  // define function f = sum over all triplets [distance from an unknown point
  // (x_0,y_0,z_0) to the plane defined by the triplet] find minimum of "f" by
  // partial derivations over x_0, y_0, and z_0 each derivation has parts
  // lineary depending on x_0, y_0, and z_0 (will fill A[deriv][3]) or to
  // nothing (will fill B[deriv]) solve A*(x_0,y_0,z_0) = B

  Acts::Vector3 vtx = Acts::Vector3::Zero();
  Acts::Vector3 vtx_prev{m_cfg.rMaxFar, m_cfg.rMaxFar, m_cfg.maxAbsZ};

  // (alpha-beta-gamma, delta), distance
  std::vector<
      std::pair<std::pair<Acts::Vector3, Acts::ActsScalar>, Acts::ActsScalar>>
      triplets_with_planes;
  triplets_with_planes.reserve(triplets.size());

  for (const auto& triplet : triplets) {
    auto abgd = makePlaneFromTriplet(triplet);
    triplets_with_planes.emplace_back(abgd, -1.);
  }

  // elements of the linear equations to solve
  Acts::SymMatrix3 A = Acts::SymMatrix3::Zero();
  Acts::Vector3 B = Acts::Vector3::Zero();
  for (const auto& triplet : triplets_with_planes) {
    const auto& abg = triplet.first.first;
    const auto& delta = triplet.first.second;

    A += 2. * (abg * abg.transpose());
    B -= 2. * delta * abg;
  }

  for (int iter = 0; iter <= m_cfg.maxIterations; iter++) {
    // new vertex position
    vtx = A.lu().solve(B);

    Acts::Vector3 vtx_diff = vtx - vtx_prev;

    if (std::sqrt(vtx_diff.dot(vtx_diff)) < m_cfg.minVtxShift) {
      // difference between the new vertex and the old vertex is not so large
      break;
    }

    if (iter != m_cfg.maxIterations) {
      // is not the last iteration
      vtx_prev = vtx;

      for (auto& triplet : triplets_with_planes) {
        const auto& abg = triplet.first.first;
        const auto& delta = triplet.first.second;
        Acts::ActsScalar distance = fabs(abg.dot(vtx) + delta);

        triplet.second = distance;
      }

      std::sort(triplets_with_planes.begin(), triplets_with_planes.end(),
                [](const auto& lhs, const auto& rhs) {
                  return lhs.second < rhs.second;
                });

      unsigned int threshold = (unsigned int)(triplets_with_planes.size() *
                                              (1. - m_cfg.removeFraction));

      for (unsigned int tr = threshold + 1; tr < triplets_with_planes.size();
           ++tr) {
        const auto& abg = triplets_with_planes.at(tr).first.first;
        const auto& delta = triplets_with_planes.at(tr).first.second;

        // remove this triplet from A and B
        A -= 2. * (abg * abg.transpose());
        B += 2. * delta * abg;
      }

      // remove all excesive triplets
      triplets_with_planes.resize(threshold);
    }
  }

  return vtx;
}

template <typename spacepoint_t>
Acts::Ray3D Acts::SeedVertexFinder<spacepoint_t>::makeRayFromTriplet(
    const Acts::SeedVertexFinder<spacepoint_t>::Triplet triplet) const {
  Acts::SymMatrix3 mat;
  mat.row(0) = Acts::Vector3(triplet.a.x(), triplet.a.y(), triplet.a.z());
  mat.row(1) = Acts::Vector3(triplet.b.x(), triplet.b.y(), triplet.b.z());
  mat.row(2) = Acts::Vector3(triplet.c.x(), triplet.c.y(), triplet.c.z());

  Acts::Vector3 mean = mat.colwise().mean();
  Acts::SymMatrix3 cov = (mat.rowwise() - mean.transpose()).transpose() *
                         (mat.rowwise() - mean.transpose()) / 3.;

  // "cov" is self-adjoint matrix
  Eigen::SelfAdjointEigenSolver<Acts::SymMatrix3> saes(cov);
  // eigenvalues are sorted in increasing order
  Acts::Vector3 eivec = saes.eigenvectors().col(2);

  return {mean, eivec};
}

template <typename spacepoint_t>
Acts::Vector3 Acts::SeedVertexFinder<spacepoint_t>::findClosestPointFromRays(
    const std::vector<Acts::SeedVertexFinder<spacepoint_t>::Triplet>& triplets)
    const {
  // define function f = sum over all triplets [distance from an unknown point
  // (x_0,y_0,z_0) to the ray defined by the triplet] find minimum of "f" by
  // partial derivations over x_0, y_0, and z_0 each derivation has parts
  // lineary depending on x_0, y_0, and z_0 (will fill A[][3]) or to nothing
  // (will fill B[]) solve A*(x_0,y_0,z_0) = B

  Acts::Vector3 vtx = Acts::Vector3::Zero();
  Acts::Vector3 vtx_prev{m_cfg.rMaxFar, m_cfg.rMaxFar, m_cfg.maxAbsZ};

  // (start_point, direction), distance
  std::vector<
      std::pair<std::pair<Acts::Vector3, Acts::Vector3>, Acts::ActsScalar>>
      triplets_with_rays;
  triplets_with_rays.reserve(triplets.size());

  for (const auto& triplet : triplets) {
    triplets_with_rays.emplace_back(
        std::make_pair(triplet.ray.origin(), triplet.ray.dir()), -1.);
  }

  // elements of the linear equations to solve
  Acts::SymMatrix3 A = Acts::SymMatrix3::Identity() * 2. * triplets.size();
  Acts::Vector3 B = Acts::Vector3::Zero();
  for (const auto& triplet : triplets_with_rays) {
    // use ray saved from earlier
    const auto& start_point = triplet.first.first;
    const auto& direction = triplet.first.second;

    A -= 2. * (direction * direction.transpose());
    B += -2. * direction * (direction.dot(start_point)) + 2. * start_point;
  }

  for (int iter = 0; iter <= m_cfg.maxIterations; iter++) {
    // new vertex position
    vtx = A.lu().solve(B);

    Acts::Vector3 vtx_diff = vtx - vtx_prev;

    if (std::sqrt(vtx_diff.dot(vtx_diff)) < m_cfg.minVtxShift) {
      // difference between the new vertex and the old vertex is not so large
      break;
    }

    if (iter != m_cfg.maxIterations) {
      // is not the last iteration
      vtx_prev = vtx;

      for (auto& triplet : triplets_with_rays) {
        const auto& start_point = triplet.first.first;
        const auto& direction = triplet.first.second;

        Acts::Vector3 vec_x_direction = (vtx - start_point).cross(direction);
        Acts::ActsScalar distance =
            std::sqrt(vec_x_direction.dot(vec_x_direction));

        triplet.second = distance;
      }

      std::sort(triplets_with_rays.begin(), triplets_with_rays.end(),
                [](const auto& lhs, const auto& rhs) {
                  return lhs.second < rhs.second;
                });

      unsigned int threshold = (unsigned int)(triplets_with_rays.size() *
                                              (1. - m_cfg.removeFraction));

      for (unsigned int tr = threshold + 1; tr < triplets_with_rays.size();
           ++tr) {
        const auto& start_point = triplets_with_rays.at(tr).first.first;
        const auto& direction = triplets_with_rays.at(tr).first.second;

        // remove this triplet from A and B
        A += 2. * (direction * direction.transpose());
        B -= -2. * direction * (direction.dot(start_point)) + 2. * start_point;
      }

      // remove all excesive triplets
      triplets_with_rays.resize(threshold);
    }
  }

  return vtx;
}
