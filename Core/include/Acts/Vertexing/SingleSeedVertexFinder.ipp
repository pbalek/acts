// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cmath>
#include <system_error>

#include <Eigen/Eigenvalues>

template <typename spacepoint_t>
Acts::SingleSeedVertexFinder<spacepoint_t>::SingleSeedVertexFinder(
    const Acts::SingleSeedVertexFinder<spacepoint_t>::Config& cfg,
    std::unique_ptr<const Logger> lgr)
    : m_cfg(cfg), m_logger(std::move(lgr)) {
  if (cfg.numPhiSlices < 3) {
    ACTS_INFO("value of numPhiSlices is "
              << cfg.numPhiSlices
              << ", which is less than 3. There will be duplicate triplets.");
  }
  if (cfg.useFracPhiSlices <= 0. || cfg.useFracPhiSlices > 1.) {
    ACTS_ERROR("value of useFracPhiSlices is "
               << cfg.useFracPhiSlices
               << ", allowed values are between 0 and 1");
  }
  if (cfg.useFracZSlices <= 0. || cfg.useFracZSlices > 1.) {
    ACTS_ERROR("value of useFracZSlices is "
               << cfg.useFracZSlices << ", allowed values are between 0 and 1");
  }
  if(cfg.minimalizeWRT != "planes" && cfg.minimalizeWRT != "rays" && cfg.minimalizeWRT != "mixed") {
    ACTS_ERROR("value of minimalizeWRT is "
               << cfg.minimalizeWRT
               << ", allowed values are \"planes\" or \"rays\" or \"mixed\".");
  }
  if(cfg.minimalizeWRT == "mixed")
  {
    if(cfg.mixedEccentricity < 0.0 || cfg.mixedEccentricity > 1.0) {
      ACTS_ERROR("value of mixedEccentricity is "
              << cfg.mixedEccentricity
              << ", allowed values are between 0 and 1.");
    }

    m_effectEccSq = cfg.mixedEccentricity*cfg.mixedEccentricity;
  }
  else {
    m_effectEccSq = (cfg.minimalizeWRT=="planes" ? 1. : 0.);
  }
  if (cfg.removeFraction < 0. || cfg.removeFraction >= 1.) {
    ACTS_ERROR("value of removeFraction is "
               << cfg.removeFraction << ", allowed values are between 0 and 1");
  }
}

template <typename spacepoint_t>
Acts::Result<Acts::Vector3>
Acts::SingleSeedVertexFinder<spacepoint_t>::findVertex(
    const std::vector<spacepoint_t>& spacepoints) const {
  ACTS_INFO("Size of spacepoints = "<<spacepoints.size()*sizeof(spacepoints[0]));
  
  // sort spacepoints to different phi and z slices
  Acts::SingleSeedVertexFinder<spacepoint_t>::SortedSpacepoints
      sortedSpacepoints = sortSpacepoints(spacepoints);

  // find triplets and fit them with plane or ray
  std::vector<Acts::SingleSeedVertexFinder<spacepoint_t>::Triplet> triplets =
      findTriplets(sortedSpacepoints);

  // if no valid triplets found
  if (triplets.empty()) {
    ACTS_INFO("triplets are empty");
    return Acts::Result<Acts::Vector3>::failure(std::error_code());
  }

  // ACTS_INFO("size of 1 Triplet "<<sizeof(Triplet)<<", size of vector7 = "<<sizeof(std::vector<double>{0.,1.,2.,3.,4.,5.,6.})+sizeof(double)*7<<" size of Ray3D "<<sizeof(Acts::Ray3D));

  // ACTS_INFO("A-Size of triplets = "<<triplets.size()<<" that is "<<triplets.size()*sizeof(triplets.at(0))<<", ActsScalar "<<sizeof(Acts::ActsScalar));

  Acts::Vector3 vtx = Acts::Vector3::Zero();

  if (m_cfg.minimalizeWRT == "planes" || m_cfg.minimalizeWRT == "rays" || m_cfg.minimalizeWRT == "mixed") {
    // find a point closest to all rays fitted through the triplets
    vtx = findClosestPoint(triplets);

    ACTS_INFO("Found mixed vertex at x = " << vtx[0]
                                     << "mm, y = " << vtx[1]
                                     << "mm, z = " << vtx[2] << "mm");

  } else {
    ACTS_ERROR("value of minimalizeWRT is "
               << m_cfg.minimalizeWRT
               << ", allowed values are \"planes\" or \"rays\" ");
  }

  // ACTS_INFO("Size of sortedSpacepoints = "<<sizeof(sortedSpacepoints.sortedSP[0][0][0])*spacepoints.size());
  // ACTS_INFO("B-Size of triplets = "<<triplets.size()<<" that is "<<triplets.size()*sizeof(triplets.at(0))*sizeof(Acts::ActsScalar)<<", ActsScalar "<<sizeof(Acts::ActsScalar));

  // FILE* file = fopen("/proc/self/status", "r");
  // char line[128];
  // while (fgets(line, 128, file) != NULL){
  //     if (strncmp(line, "VmSize:", 7) == 0){
  //         ACTS_INFO("In the end of SingleSeedVertexFinder: VmSize="<<line);
  //     }
  //     if (strncmp(line, "VmRSS:", 6) == 0){
  //       ACTS_INFO("In the end of SingleSeedVertexFinder: VmRSS="<<line);
  //     }
  // }
  // fclose(file);

  return Acts::Result<Acts::Vector3>::success(vtx);
}

template <typename spacepoint_t>
typename Acts::SingleSeedVertexFinder<spacepoint_t>::SortedSpacepoints
Acts::SingleSeedVertexFinder<spacepoint_t>::sortSpacepoints(
    const std::vector<spacepoint_t>& spacepoints) const {
  Acts::SingleSeedVertexFinder<spacepoint_t>::SortedSpacepoints
      sortedSpacepoints(m_cfg.numPhiSlices, m_cfg.numZSlices);

  for (const auto& sp : spacepoints) {
    // phi will be saved for later
    Acts::ActsScalar phi = detail::radian_pos(std::atan2(sp.y(), sp.x()));
    std::uint32_t phislice =
        (std::uint32_t)(phi / (2 * M_PI) * m_cfg.numPhiSlices);
    if (phislice >= m_cfg.numPhiSlices) {
      phislice = 0;
    }

    if (std::abs(sp.z()) >= m_cfg.maxAbsZ) {
      continue;
    }
    std::uint32_t zslice = (std::uint32_t)(
        (sp.z() + m_cfg.maxAbsZ) / (2 * m_cfg.maxAbsZ) * m_cfg.numZSlices);

    // input spacepoint is sorted into one of the subsets
    if (sp.r() < m_cfg.rMinMiddle) {
      if (m_cfg.rMinNear < sp.r() && sp.r() < m_cfg.rMaxNear) {
        if (std::fmod(m_cfg.useFracPhiSlices * phislice, 1.0) >=
            m_cfg.useFracPhiSlices) {
          continue;
        }
        sortedSpacepoints.addSP(0, phislice, zslice)
            .emplace_back((spacepoint_t const*)&sp, phi);
      }
    } else if (sp.r() < m_cfg.rMinFar) {
      if (sp.r() < m_cfg.rMaxMiddle) {
        if (std::fmod(m_cfg.useFracZSlices * zslice, 1.0) >=
            m_cfg.useFracZSlices) {
          continue;
        }
        sortedSpacepoints.addSP(1, phislice, zslice)
            .emplace_back((spacepoint_t const*)&sp, phi);
      }
    } else if (sp.r() < m_cfg.rMaxFar) {
      sortedSpacepoints.addSP(2, phislice, zslice)
          .emplace_back((spacepoint_t const*)&sp, phi);
    }
  }

  std::uint32_t counter[3]={0,0,0};
  for(std::uint32_t la=0;la<3;++la)
  {
    for(std::uint32_t phi=0;phi<m_cfg.numPhiSlices;++phi)
    {
      for(std::uint32_t z=0;z<m_cfg.numZSlices;++z)
      {
        counter[la]+=sortedSpacepoints.getSP(la,phi,z).size();
      }
    }
  }
  ACTS_INFO("have sorted spacepoints: "<<counter[0]<<", "<<counter[1]<<", "<<counter[2]);

  return sortedSpacepoints;
}

template <typename spacepoint_t>
std::vector<typename Acts::SingleSeedVertexFinder<spacepoint_t>::Triplet>
Acts::SingleSeedVertexFinder<spacepoint_t>::findTriplets(
    const Acts::SingleSeedVertexFinder<spacepoint_t>::SortedSpacepoints&
        sortedSpacepoints) const {
  std::vector<Acts::SingleSeedVertexFinder<spacepoint_t>::Triplet> triplets;

  std::uint32_t phiStep =
      (std::uint32_t)(m_cfg.maxPhideviation / (2 * M_PI / m_cfg.numPhiSlices)) +
      1;

  std::uint32_t considered_triplets = 0;
  std::uint32_t rej0=0, rej1 = 0, rej2 = 0, rej3 = 0, rej4 = 0;
  std::uint32_t good_triplets = 0;

  // calculate limits for middle spacepoints
  Acts::Vector2 vecA{-m_cfg.maxAbsZ + m_cfg.maxZPosition, m_cfg.rMinFar};
  vecA /= 2.;
  Acts::Vector2 vecB = {vecA[1], -vecA[0]};
  vecB /= std::tan(m_cfg.maxXYZdeviation);
  Acts::Vector2 posR = Acts::Vector2(-m_cfg.maxZPosition, 0.) + vecA + vecB;
  Acts::ActsScalar R = vecA.norm() / std::sin(m_cfg.maxXYZdeviation);
  Acts::ActsScalar constB = -2. * posR[0];
  Acts::ActsScalar constC =
      posR[0] * posR[0] +
      (posR[1] - m_cfg.rMaxNear) * (posR[1] - m_cfg.rMaxNear) - R * R;
  Acts::ActsScalar maxZMiddle =
      -1. * (-constB - sqrt(constB * constB - 4. * constC)) / 2.;
  if (maxZMiddle <= 0) {
    ACTS_WARNING(
        "maximum position of middle spacepoints is not positive, maxZMiddle = "
        << maxZMiddle << ", check your config; setting maxZMiddle to "
        << m_cfg.maxAbsZ);
    maxZMiddle = m_cfg.maxAbsZ;
  }

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
  std::uint32_t limitMiddleSliceFrom =
      (std::uint32_t)((-maxZMiddle + m_cfg.maxAbsZ) / zBinLength);
  std::uint32_t limitMiddleSliceTo =
      (std::uint32_t)((maxZMiddle + m_cfg.maxAbsZ) / zBinLength + 1);
  std::uint32_t limitAbsZSliceFrom = (std::uint32_t)(
      (-m_cfg.maxZPosition + m_cfg.maxAbsZ) / zBinLength + 0.01);
  std::uint32_t limitAbsZSliceTo =
      (std::uint32_t)((m_cfg.maxZPosition + m_cfg.maxAbsZ) / zBinLength + 1.01);

  for (std::uint32_t middleZ = limitMiddleSliceFrom;
       middleZ < limitMiddleSliceTo; ++middleZ) {
    // skip slices that are empty anyway
    if (std::fmod(m_cfg.useFracZSlices * middleZ, 1.0) >=
        m_cfg.useFracZSlices) {
      continue;
    }

    // calculate limits for near spacepoints, assuming the middle spacepoints
    // are within some boundaries
    bool isLessFrom = (middleZ <= limitAbsZSliceFrom);
    Acts::ActsScalar deltaZfrom =
        (middleZ - limitAbsZSliceFrom - 1) * zBinLength;
    Acts::ActsScalar angleZfrom =
        std::atan2(rMiddle[isLessFrom], deltaZfrom) + m_cfg.maxXYZdeviation;
    std::uint32_t nearZFrom = 0;
    if (angleZfrom < M_PI) {
      Acts::ActsScalar new_deltaZfrom =
          rMiddle[isLessFrom] / std::tan(angleZfrom) / zBinLength;
      nearZFrom = (std::uint32_t)std::max(
          new_deltaZfrom * rNearRatio[isLessFrom] + limitAbsZSliceFrom, 0.);
    }

    bool isLessTo = (middleZ < limitAbsZSliceTo);
    Acts::ActsScalar deltaZto = (middleZ - limitAbsZSliceTo + 1) * zBinLength;
    Acts::ActsScalar angleZto =
        std::atan2(rMiddle[!isLessTo], deltaZto) - m_cfg.maxXYZdeviation;
    std::uint32_t nearZTo = m_cfg.numZSlices;
    if (angleZto > 0) {
      Acts::ActsScalar new_deltaZto =
          rMiddle[!isLessTo] / std::tan(angleZto) / zBinLength;
      nearZTo = (std::uint32_t)std::max(
          new_deltaZto * rNearRatio[!isLessTo] + limitAbsZSliceTo, 0.);
      if (nearZTo > m_cfg.numZSlices) {
        nearZTo = m_cfg.numZSlices;
      }
    }

    for (std::uint32_t nearZ = nearZFrom; nearZ < nearZTo; ++nearZ) {
      // calculate limits for far spacepoits, assuming middle and near
      // spacepoits are within some boundaries
      bool isMiddleLess = (middleZ <= nearZ);

      Acts::ActsScalar delta2Zfrom = (middleZ - nearZ - 1) * zBinLength;
      Acts::ActsScalar angle2Zfrom =
          std::atan2(rFarDelta[isMiddleLess], delta2Zfrom) +
          m_cfg.maxXYZdeviation;
      std::uint32_t farZFrom = 0;
      if (angle2Zfrom < M_PI) {
        farZFrom = (std::uint32_t)std::max(
            (rFarDelta[isMiddleLess] / std::tan(angle2Zfrom) / zBinLength) +
                middleZ,
            0.);
        if (farZFrom >= m_cfg.numZSlices) {
          continue;
        }
      }

      isMiddleLess = (middleZ < nearZ);
      Acts::ActsScalar delta2Zto = (middleZ - nearZ + 1) * zBinLength;
      Acts::ActsScalar angle2Zto =
          std::atan2(rFarDelta[!isMiddleLess], delta2Zto) -
          m_cfg.maxXYZdeviation;
      std::uint32_t farZTo = m_cfg.numZSlices;
      if (angle2Zto > 0) {
        farZTo = (std::uint32_t)std::max(
            (rFarDelta[!isMiddleLess] / std::tan(angle2Zto) / zBinLength) +
                middleZ + 1,
            0.);
        if (farZTo > m_cfg.numZSlices) {
          farZTo = m_cfg.numZSlices;
        } else if (farZTo == 0) {
          continue;
        }
      }

      for (std::uint32_t farZ = farZFrom; farZ < farZTo; farZ++) {
        // loop over near phi slices
        for (std::uint32_t nearPhi = 0; nearPhi < m_cfg.numPhiSlices;
             ++nearPhi) {
          // skip slices that are empty anyway
          if (std::fmod(m_cfg.useFracPhiSlices * nearPhi, 1.0) >=
              m_cfg.useFracPhiSlices) {
            continue;
          }

          // loop over some middle phi slices
          for (std::uint32_t middlePhi_h =
                   m_cfg.numPhiSlices + nearPhi - phiStep;
               middlePhi_h <= m_cfg.numPhiSlices + nearPhi + phiStep;
               ++middlePhi_h) {
            std::uint32_t middlePhi = middlePhi_h % m_cfg.numPhiSlices;

            // loop over some far phi slices
            for (std::uint32_t farPhi_h =
                     m_cfg.numPhiSlices + middlePhi - phiStep;
                 farPhi_h <= m_cfg.numPhiSlices + middlePhi + phiStep;
                 ++farPhi_h) {
              std::uint32_t farPhi = farPhi_h % m_cfg.numPhiSlices;

              // for all near spacepoints in this slice
              for (const auto& nearSP :
                   sortedSpacepoints.getSP(0, nearPhi, nearZ)) {
                Acts::ActsScalar phiA = nearSP.second;

                // for all middle spacepoints in this slice
                for (const auto& middleSP :
                     sortedSpacepoints.getSP(1, middlePhi, middleZ)) {
                  Acts::ActsScalar phiB = middleSP.second;
                  Acts::ActsScalar deltaPhiAB =
                      detail::difference_periodic(phiA, phiB, 2 * M_PI);
                  if (std::abs(deltaPhiAB) > m_cfg.maxPhideviation) {
                    continue;
                  }

                  // for all far spacepoints in this slice
                  for (const auto& farSP :
                       sortedSpacepoints.getSP(2, farPhi, farZ)) {
                    Acts::ActsScalar phiC = farSP.second;
                    Acts::ActsScalar deltaPhiBC =
                        detail::difference_periodic(phiB, phiC, 2 * M_PI);
                    if (std::abs(deltaPhiBC) > m_cfg.maxPhideviation) {
                      continue;
                    }

                    auto fittedTriplet = tripletValidationAndFit(
                        *nearSP.first, *middleSP.first, *farSP.first);
                    
                    ++considered_triplets;

                    if (fittedTriplet.getDistance()>-1.5) {
                      triplets.push_back(fittedTriplet);
                      ++good_triplets;
                    }
                    else
                    {
                      ++rej0;
                      if(fittedTriplet.getDistance()>-2.15) ++rej1;
                      if(fittedTriplet.getDistance()>-2.25) ++rej2;
                      if(fittedTriplet.getDistance()>-2.35) ++rej3;
                      if(fittedTriplet.getDistance()>-2.45) ++rej4;
                    }
                  }  // loop over far spacepoints
                }    // loop over middle spacepoints
              }      // loop over near spacepoints
            }        // loop over far phi slices
          }          // loop over middle phi slices
        }            // loop over near phi slices
      }              // loop over far Z slices
    }                // loop over near Z slices
  }                  // loop over middle Z slices

  ACTS_INFO("considered "<<considered_triplets<<", good "<<good_triplets<<", i.e. bad "<<considered_triplets-good_triplets);
  ACTS_INFO("rejected  total - "<<rej0<<"; rejected by step 1 - "<<rej1<<", by step 2 - "<<rej2<<", by step 3 - "<<rej3<<", by step 4 - "<<rej4);
  return triplets;
}

template <typename spacepoint_t>
typename Acts::SingleSeedVertexFinder<spacepoint_t>::Triplet Acts::SingleSeedVertexFinder<spacepoint_t>::tripletValidationAndFit(
    const spacepoint_t& a, const spacepoint_t& b, const spacepoint_t& c) const {
  // slope for near+middle spacepoints
  Acts::ActsScalar alpha1 =
      std::atan2(a.y() - b.y(), a.x() - b.x());
  // slope for middle+far spacepoints
  Acts::ActsScalar alpha2 =
      std::atan2(b.y() - c.y(), b.x() - c.x());
  // these two slopes shouldn't be too different
  Acts::ActsScalar deltaAlpha =
      detail::difference_periodic(alpha1, alpha2, 2 * M_PI);
  if (std::abs(deltaAlpha) > m_cfg.maxXYdeviation) {
    return Triplet(-2.1, a,b,c);
  }

  // near-middle ray
  Acts::Vector3 ab{a.x() - b.x(), a.y() - b.y(),
                   a.z() - b.z()};
  // middle-far ray
  Acts::Vector3 bc{b.x() - c.x(), b.y() - c.y(),
                   b.z() - c.z()};
  // dot product of these two
  Acts::ActsScalar cosTheta = (ab.dot(bc)) / (ab.norm() * bc.norm());
  Acts::ActsScalar theta = std::acos(cosTheta);
  if (theta > m_cfg.maxXYZdeviation) {
    return Triplet(-2.2, a,b,c);
  }

  // reject the ray if it doesn't come close to the z-axis
  auto triplet = fitTriplet(a, b, c);
  const Acts::Vector3& startPoint = triplet.getStartPoint();
  const Acts::Vector3& direction = triplet.getDirection();
  // ACTS_INFO("fitted triplet: "<<startPoint[0]<<" "<<startPoint[1]<<" "<<startPoint[2]<<" .. "<<direction[0]<<" "<<direction[1]<<" "<<direction[2]);
  // norm to z-axis and to the ray
  Acts::Vector3 norm{-1. * direction[1], 1. * direction[0], 0};
  Acts::ActsScalar norm_size = norm.norm();

  Acts::ActsScalar tanTheta = norm_size / direction[2];
  if (std::abs(tanTheta) < std::tan(m_cfg.minTheta)) {
    return Triplet(-2.3,a,b,c);
  }

  // nearest distance from the ray to z-axis
  Acts::ActsScalar dist = std::abs(startPoint.dot(norm)) / norm_size;
  if (dist > m_cfg.maxRPosition) {
    return Triplet(-2.4, a,b,c);
  }

  // z coordinate of the nearest distance from the ray to z-axis
  Acts::ActsScalar zDist =
      direction.cross(norm).dot(startPoint) / (norm_size * norm_size);
  if (std::abs(zDist) > m_cfg.maxZPosition) {
    return Triplet(-2.5, a,b,c);
  }

  return triplet;
}

template <typename spacepoint_t>
typename Acts::SingleSeedVertexFinder<spacepoint_t>::Triplet
Acts::SingleSeedVertexFinder<spacepoint_t>::fitTriplet(
    const spacepoint_t& spA, const spacepoint_t& spB, const spacepoint_t& spC) const {
  // always fit ray
  auto ray = makeRayFromTriplet(spA,spB,spC);

  std::pair<Acts::Vector3, Acts::ActsScalar> plane{{0.,0.,1.},0.};
  if(m_effectEccSq!=0.)
  {
    // fit plane only if minimalizeWRT is not "rays"
    plane = makePlaneFromTriplet(spA,spB,spC);
  }

  return Triplet(plane.first, plane.second, ray.first, ray.second, spA,spB,spC);
}

template <typename spacepoint_t>
std::pair<Acts::Vector3, Acts::ActsScalar>
Acts::SingleSeedVertexFinder<spacepoint_t>::makePlaneFromTriplet(
    const spacepoint_t& spA, const spacepoint_t& spB, const spacepoint_t& spC) const {
  Acts::Vector3 a{spA.x(), spA.y(), spA.z()};
  Acts::Vector3 b{spB.x(), spB.y(), spB.z()};
  Acts::Vector3 c{spC.x(), spC.y(), spC.z()};

  Acts::Vector3 ba = b - a, ca = c - a;

  // vector (alpha,beta,gamma) normalized to unity for convenience
  Acts::Vector3 abg = ba.cross(ca).normalized();
  Acts::ActsScalar delta = -1. * abg.dot(a);

  // plane (alpha*x + beta*y + gamma*z + delta = 0)
  return {abg, delta};
}

template <typename spacepoint_t>
std::pair<Acts::Vector3, Acts::Vector3> Acts::SingleSeedVertexFinder<spacepoint_t>::makeRayFromTriplet(
    const spacepoint_t& a, const spacepoint_t& b, const spacepoint_t& c) const {
  Acts::SymMatrix3 mat;
  mat.row(0) = Acts::Vector3(a.x(), a.y(), a.z());
  mat.row(1) = Acts::Vector3(b.x(), b.y(), b.z());
  mat.row(2) = Acts::Vector3(c.x(), c.y(), c.z());

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
Acts::Vector3
Acts::SingleSeedVertexFinder<spacepoint_t>::findClosestPoint(
    std::vector<typename Acts::SingleSeedVertexFinder<spacepoint_t>::Triplet>& allTriples) const {
  // 1. define function f = sum over all triplets [distance from an unknown
  // point
  //    (x_0,y_0,z_0) to the plane defined by the triplet]
  // 2. find minimum of "f" by partial derivations over x_0, y_0, and z_0
  // 3. each derivation has parts linearly depending on x_0, y_0, and z_0
  //    (will fill A[deriv][3]) or to nothing (will fill B[deriv])
  // 4. solve A*(x_0,y_0,z_0) = B

  Acts::Vector3 vtx = Acts::Vector3::Zero();
  Acts::Vector3 vtxPrev{m_cfg.rMaxFar, m_cfg.rMaxFar, m_cfg.maxAbsZ};


  // ACTS_INFO("A-size of Acts::Vector3 "<<sizeof(Acts::Vector3)<<", "<<sizeof(std::pair<Acts::Vector3, Acts::ActsScalar>)<<"; tripletsWithRays pairs "<<sizeof(std::pair<std::pair<Acts::Vector3, Acts::ActsScalar>, Acts::ActsScalar>)<<", size "<<allTriples.size());

  ACTS_INFO("A-m_effectEccSq "<<m_effectEccSq);

  // elements of the linear equations to solve
  Acts::SymMatrix3 A = Acts::SymMatrix3::Zero() + (1.-m_effectEccSq)*Acts::SymMatrix3::Identity() * 2. * allTriples.size();
  Acts::Vector3 B = Acts::Vector3::Zero();
  for (const auto& triplet : allTriples) {
    const Acts::Vector3& abg = triplet.getPlaneABG();
    const Acts::ActsScalar& delta = triplet.getPlaneDelta();
    const Acts::Vector3& startPoint=triplet.getStartPoint();
    const Acts::Vector3& direction=triplet.getDirection();

    A += m_effectEccSq * 2. * (abg * abg.transpose()) - (1-m_effectEccSq) * 2. * (direction * direction.transpose());
    B -= m_effectEccSq * 2. * delta * abg + (1-m_effectEccSq) * (2. * direction * (direction.dot(startPoint)) - 2. * startPoint);
  }

  for (std::uint32_t iter = 0; iter <= m_cfg.maxIterations; iter++) {
    // new vertex position
    vtx = A.lu().solve(B);

    Acts::Vector3 vtxDiff = vtx - vtxPrev;

    ACTS_INFO(iter<<": vtx = "<<vtx[0]<<", "<<vtx[1]<<", "<<vtx[2]<<"; vtxDiff = "<<vtxDiff[0]<<", "<<vtxDiff[1]<<", "<<vtxDiff[2]<<"; there are "<<allTriples.size()<<" triplets");

    if (vtxDiff.norm() < m_cfg.minVtxShift) {
      // difference between the new vertex and the old vertex is not so large
      break;
    }

    if (iter != m_cfg.maxIterations) {
      // is not the last iteration
      vtxPrev = vtx;

      for (auto& triplet : allTriples) {
        const Acts::Vector3& abg = triplet.getPlaneABG();
        const Acts::ActsScalar& delta = triplet.getPlaneDelta();
        const Acts::Vector3& startPoint=triplet.getStartPoint();
        const Acts::Vector3& direction=triplet.getDirection();

        const Acts::ActsScalar distanceSq = m_effectEccSq * std::abs(abg.dot(vtx) + delta) * std::abs(abg.dot(vtx) + delta) + (1.-m_effectEccSq) * (vtx - startPoint).cross(direction).squaredNorm();
        triplet.setDistance(distanceSq);
      }

      std::sort(allTriples.begin(), allTriples.end(),
                [](const auto& lhs, const auto& rhs) {
                  return lhs.getDistance() < rhs.getDistance();
                });

      std::uint32_t threshold = (std::uint32_t)(allTriples.size() *
                                                (1. - m_cfg.removeFraction));

      for (std::uint32_t tr = threshold + 1; tr < allTriples.size();
           ++tr) {
        const Acts::Vector3& abg = allTriples[tr].getPlaneABG();
        const Acts::ActsScalar& delta = allTriples[tr].getPlaneDelta();
        const Acts::Vector3& startPoint = allTriples[tr].getStartPoint();
        const Acts::Vector3& direction  = allTriples[tr].getDirection();

        // remove this triplet from A and B
        A -= m_effectEccSq * 2. * (abg * abg.transpose()) + (1-m_effectEccSq) * (Acts::SymMatrix3::Identity() * 2. -  2. * (direction * direction.transpose()));
        B += m_effectEccSq * 2. * delta * abg + (1-m_effectEccSq) * (2. * direction * (direction.dot(startPoint)) - 2. * startPoint);
      }

      // remove all excessive triplets
      if(threshold + 1 < allTriples.size()) {
        ACTS_INFO("removed "<<allTriples.size()-threshold-1<<" triplets with removeFraction");
        allTriples.erase(allTriples.begin()+threshold+1,allTriples.end());
      }

      // std::uint32_t remove_counter=0;
      // std::uint32_t good_threshold = (std::uint32_t)(tripletsWithRays.size()*0.40);
      // std::uint32_t bad_threshold = (std::uint32_t)(tripletsWithRays.size()*0.80);
      // for (std::uint32_t tr = 0; tr < good_threshold; ++tr) {

      //   auto& good_triplet = tripletsWithRays[tr];

      //   for (std::uint32_t tr2 = std::max(tr+1,bad_threshold); tr2 < tripletsWithRays.size(); ++tr2) {
      //     auto& bad_triplet = tripletsWithRays[tr2];

      //     if(good_triplet.nearSP==bad_triplet.nearSP || 
      //        good_triplet.middleSP==bad_triplet.middleSP || 
      //        good_triplet.farSP==bad_triplet.farSP) {
      //       //ACTS_INFO("good triplet tr "<<tr<<" has overlap with bad triplet tr2 "<<tr2);

      //       const Acts::Vector3& abg = tripletsWithPlanes[tr2].getPlaneABG();
      //       const Acts::ActsScalar& delta = tripletsWithPlanes[tr2].getPlaneDelta();

      //       // remove this triplet from A and B
      //       A -= 2. * (abg * abg.transpose());
      //       B += 2. * delta * abg;

      //       tripletsWithPlanes.erase(tripletsWithPlanes.begin()+tr2,tripletsWithPlanes.begin()+tr2+1);

      //       --tr2;
      //       ++remove_counter;
      //     }
      //   }
      // }
      // ACTS_INFO("removed "<<remove_counter<<" triplets with good-bad triplets / 1");

      // remove_counter=0;
      // good_threshold = (std::uint32_t)(tripletsWithPlanes.size()*0.90);
      // bad_threshold = (std::uint32_t)(tripletsWithPlanes.size()*0.60);
      // for (std::uint32_t tr = 0; tr < good_threshold; ++tr) {

      //   auto& good_triplet = tripletsWithPlanes[tr];

      //   for (std::uint32_t tr2 = std::max(tr+1,bad_threshold); tr2 < tripletsWithPlanes.size(); ++tr2) {
      //     auto& bad_triplet = tripletsWithPlanes[tr2];

      //     if((good_triplet.nearSP==bad_triplet.nearSP && 
      //         (good_triplet.middleSP==bad_triplet.middleSP || good_triplet.farSP==bad_triplet.farSP)) ||
      //         (good_triplet.middleSP==bad_triplet.middleSP && good_triplet.farSP==bad_triplet.farSP)) {
      //       //ACTS_INFO("good triplet tr "<<tr<<" has overlap with bad triplet tr2 "<<tr2);

      //       const Acts::Vector3& abg = tripletsWithPlanes[tr2].getPlaneABG();
      //       const Acts::ActsScalar& delta = tripletsWithPlanes[tr2].getPlaneDelta();

      //       // remove this triplet from A and B
      //       A -= 2. * (abg * abg.transpose());
      //       B += 2. * delta * abg;

      //       tripletsWithPlanes.erase(tripletsWithPlanes.begin()+tr2,tripletsWithPlanes.begin()+tr2+1);

      //       --tr2;
      //       ++remove_counter;
      //     }
      //   }
      // }
      // ACTS_INFO("removed "<<remove_counter<<" triplets with good-bad triplets / 2");


    }
  }

  ACTS_INFO("B-size of Acts::Vector3 "<<sizeof(Acts::Vector3)<<", "<<sizeof(std::pair<Acts::Vector3, Acts::ActsScalar>)<<"; allTriples pairs "<<sizeof(std::pair<std::pair<Acts::Vector3, Acts::ActsScalar>, Acts::ActsScalar>)<<", size "<<allTriples.size());

  return vtx;
}

////////////////////////////////

/*
template <typename spacepoint_t>
Acts::Vector3
Acts::SingleSeedVertexFinder<spacepoint_t>::findClosestPointFromPlanes(
    std::vector<typename Acts::SingleSeedVertexFinder<spacepoint_t>::Triplet>& tripletsWithPlanes) const {
  // 1. define function f = sum over all triplets [distance from an unknown
  // point
  //    (x_0,y_0,z_0) to the plane defined by the triplet]
  // 2. find minimum of "f" by partial derivations over x_0, y_0, and z_0
  // 3. each derivation has parts linearly depending on x_0, y_0, and z_0
  //    (will fill A[deriv][3]) or to nothing (will fill B[deriv])
  // 4. solve A*(x_0,y_0,z_0) = B

  Acts::Vector3 vtx = Acts::Vector3::Zero();
  Acts::Vector3 vtxPrev{m_cfg.rMaxFar, m_cfg.rMaxFar, m_cfg.maxAbsZ};


  ACTS_INFO("A-size of Acts::Vector3 "<<sizeof(Acts::Vector3)<<", "<<sizeof(std::pair<Acts::Vector3, Acts::ActsScalar>)<<"; tripletsWithRays pairs "<<sizeof(std::pair<std::pair<Acts::Vector3, Acts::ActsScalar>, Acts::ActsScalar>)<<", size "<<tripletsWithPlanes.size());


  // elements of the linear equations to solve
  Acts::SymMatrix3 A = Acts::SymMatrix3::Zero();
  Acts::Vector3 B = Acts::Vector3::Zero();
  for (const auto& triplet : tripletsWithPlanes) {
    const Acts::Vector3& abg = triplet.getPlaneABG();
    const Acts::ActsScalar& delta = triplet.getPlaneDelta();

    A += 2. * (abg * abg.transpose());
    B -= 2. * delta * abg;
  }

  for (std::uint32_t iter = 0; iter <= m_cfg.maxIterations; iter++) {
    // new vertex position
    vtx = A.lu().solve(B);

    Acts::Vector3 vtxDiff = vtx - vtxPrev;

    ACTS_INFO(iter<<": vtx = "<<vtx[0]<<", "<<vtx[1]<<", "<<vtx[2]<<"; vtxDiff = "<<vtxDiff[0]<<", "<<vtxDiff[1]<<", "<<vtxDiff[2]<<"; there are "<<tripletsWithPlanes.size()<<" triplets");

    if (vtxDiff.norm() < m_cfg.minVtxShift) {
      // difference between the new vertex and the old vertex is not so large
      break;
    }

    if (iter != m_cfg.maxIterations) {
      // is not the last iteration
      vtxPrev = vtx;

      for (auto& triplet : tripletsWithPlanes) {
        const Acts::Vector3& abg = triplet.getPlaneABG();
        const Acts::ActsScalar& delta = triplet.getPlaneDelta();
        const Acts::ActsScalar distance = std::abs(abg.dot(vtx) + delta);
        triplet.setDistance(distance);
      }

      std::sort(tripletsWithPlanes.begin(), tripletsWithPlanes.end(),
                [](const auto& lhs, const auto& rhs) {
                  return lhs.getDistance() < rhs.getDistance();
                });

      std::uint32_t threshold = (std::uint32_t)(tripletsWithPlanes.size() *
                                                (1. - m_cfg.removeFraction));

      for (std::uint32_t tr = threshold + 1; tr < tripletsWithPlanes.size();
           ++tr) {
        const Acts::Vector3& abg = tripletsWithPlanes[tr].getPlaneABG();
        const Acts::ActsScalar& delta = tripletsWithPlanes[tr].getPlaneDelta();

        // remove this triplet from A and B
        A -= 2. * (abg * abg.transpose());
        B += 2. * delta * abg;
      }

      // remove all excessive triplets
      if(threshold + 1 < tripletsWithPlanes.size()) {
        ACTS_INFO("removed "<<tripletsWithPlanes.size()-threshold-1<<" triplets with removeFraction");
        tripletsWithPlanes.erase(tripletsWithPlanes.begin()+threshold+1,tripletsWithPlanes.end());
      }
*/
/*
      std::uint32_t remove_counter=0;
      std::uint32_t good_threshold = (std::uint32_t)(tripletsWithPlanes.size()*0.40);
      std::uint32_t bad_threshold = (std::uint32_t)(tripletsWithPlanes.size()*0.80);
      for (std::uint32_t tr = 0; tr < good_threshold; ++tr) {

        auto& good_triplet = tripletsWithPlanes[tr];

        for (std::uint32_t tr2 = std::max(tr+1,bad_threshold); tr2 < tripletsWithPlanes.size(); ++tr2) {
          auto& bad_triplet = tripletsWithPlanes[tr2];

          if(good_triplet.nearSP==bad_triplet.nearSP || 
             good_triplet.middleSP==bad_triplet.middleSP || 
             good_triplet.farSP==bad_triplet.farSP) {
            //ACTS_INFO("good triplet tr "<<tr<<" has overlap with bad triplet tr2 "<<tr2);

            const Acts::Vector3& abg = tripletsWithPlanes[tr2].getPlaneABG();
            const Acts::ActsScalar& delta = tripletsWithPlanes[tr2].getPlaneDelta();

            // remove this triplet from A and B
            A -= 2. * (abg * abg.transpose());
            B += 2. * delta * abg;

            tripletsWithPlanes.erase(tripletsWithPlanes.begin()+tr2,tripletsWithPlanes.begin()+tr2+1);

            --tr2;
            ++remove_counter;
          }
        }
      }
      ACTS_INFO("removed "<<remove_counter<<" triplets with good-bad triplets / 1");

      remove_counter=0;
      good_threshold = (std::uint32_t)(tripletsWithPlanes.size()*0.90);
      bad_threshold = (std::uint32_t)(tripletsWithPlanes.size()*0.60);
      for (std::uint32_t tr = 0; tr < good_threshold; ++tr) {

        auto& good_triplet = tripletsWithPlanes[tr];

        for (std::uint32_t tr2 = std::max(tr+1,bad_threshold); tr2 < tripletsWithPlanes.size(); ++tr2) {
          auto& bad_triplet = tripletsWithPlanes[tr2];

          if((good_triplet.nearSP==bad_triplet.nearSP && 
              (good_triplet.middleSP==bad_triplet.middleSP || good_triplet.farSP==bad_triplet.farSP)) ||
              (good_triplet.middleSP==bad_triplet.middleSP && good_triplet.farSP==bad_triplet.farSP)) {
            //ACTS_INFO("good triplet tr "<<tr<<" has overlap with bad triplet tr2 "<<tr2);

            const Acts::Vector3& abg = tripletsWithPlanes[tr2].getPlaneABG();
            const Acts::ActsScalar& delta = tripletsWithPlanes[tr2].getPlaneDelta();

            // remove this triplet from A and B
            A -= 2. * (abg * abg.transpose());
            B += 2. * delta * abg;

            tripletsWithPlanes.erase(tripletsWithPlanes.begin()+tr2,tripletsWithPlanes.begin()+tr2+1);

            --tr2;
            ++remove_counter;
          }
        }
      }
      ACTS_INFO("removed "<<remove_counter<<" triplets with good-bad triplets / 2");
*/
/*
    }
  }

  ACTS_INFO("B-size of Acts::Vector3 "<<sizeof(Acts::Vector3)<<", "<<sizeof(std::pair<Acts::Vector3, Acts::ActsScalar>)<<"; tripletsWithRays pairs "<<sizeof(std::pair<std::pair<Acts::Vector3, Acts::ActsScalar>, Acts::ActsScalar>)<<", size "<<tripletsWithPlanes.size());

  return vtx;
}



template <typename spacepoint_t>
Acts::Vector3
Acts::SingleSeedVertexFinder<spacepoint_t>::findClosestPointFromRays(
    std::vector<typename Acts::SingleSeedVertexFinder<spacepoint_t>::Triplet>&
        tripletsWithRays) const {
  // 1. define function f = sum over all triplets [distance from an unknown
  // point
  //    (x_0,y_0,z_0) to the ray defined by the triplet]
  // 2. find minimum of "f" by partial derivations over x_0, y_0, and z_0
  // 3. each derivation has parts linearly depending on x_0, y_0, and z_0
  //    (will fill A[][3]) or to nothing (will fill B[])
  // 4. solve A*(x_0,y_0,z_0) = B


  Acts::Vector3 vtx = Acts::Vector3::Zero();
  Acts::Vector3 vtxPrev{m_cfg.rMaxFar, m_cfg.rMaxFar, m_cfg.maxAbsZ};

  // elements of the linear equations to solve
  Acts::SymMatrix3 A = Acts::SymMatrix3::Identity() * 2. * tripletsWithRays.size();
  Acts::Vector3 B = Acts::Vector3::Zero();
  for (const auto& triplet : tripletsWithRays) {
    // use ray saved from earlier
    const Acts::Vector3& startPoint=triplet.getStartPoint();
    const Acts::Vector3& direction=triplet.getDirection();

    A -= 2. * (direction * direction.transpose());
    B += -2. * direction * (direction.dot(startPoint)) + 2. * startPoint;
  }

  for (std::uint32_t iter = 0; iter <= m_cfg.maxIterations; iter++) {
    // new vertex position
    vtx = A.lu().solve(B);

    Acts::Vector3 vtxDiff = vtx - vtxPrev;

    ACTS_INFO(iter<<": RAY vtx = "<<vtx[0]<<", "<<vtx[1]<<", "<<vtx[2]<<"; vtxDiff = "<<vtxDiff[0]<<", "<<vtxDiff[1]<<", "<<vtxDiff[2]<<"; there are "<<tripletsWithRays.size()<<" triplets");

    if (vtxDiff.norm() < m_cfg.minVtxShift) {
      // difference between the new vertex and the old vertex is not so large
      break;
    }

    if (iter != m_cfg.maxIterations) {
      // is not the last iteration
      vtxPrev = vtx;

      for (auto& triplet : tripletsWithRays) {
        const Acts::Vector3& startPoint=triplet.getStartPoint();
        const Acts::Vector3& direction=triplet.getDirection();
        const Acts::ActsScalar distance = (vtx - startPoint).cross(direction).norm();
        triplet.setDistance(distance);
      }

      std::sort(tripletsWithRays.begin(), tripletsWithRays.end(),
                [](const auto& lhs, const auto& rhs) {
                  return lhs.getDistance() < rhs.getDistance();
                });

      std::uint32_t threshold = (std::uint32_t)(tripletsWithRays.size() *
                                                (1. - m_cfg.removeFraction));

      for (std::uint32_t tr = threshold + 1; tr < tripletsWithRays.size();
           ++tr) {
        const Acts::Vector3& startPoint = tripletsWithRays[tr].getStartPoint();
        const Acts::Vector3& direction  = tripletsWithRays[tr].getDirection();

        // remove this triplet from A and B
        A -= Acts::SymMatrix3::Identity() * 2.;
        A += 2. * (direction * direction.transpose());
        B -= -2. * direction * (direction.dot(startPoint)) + 2. * startPoint;
      }

      // remove all excessive triplets
      if(threshold + 1 < tripletsWithRays.size()) {
        ACTS_INFO("removed "<<tripletsWithRays.size()-threshold-1<<" triplets with removeFraction");
        tripletsWithRays.erase(tripletsWithRays.begin()+threshold+1,tripletsWithRays.end());
      }
    }
  }

  return vtx;
}
*/