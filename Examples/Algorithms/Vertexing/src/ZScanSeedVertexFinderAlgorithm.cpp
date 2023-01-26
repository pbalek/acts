// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Vertexing/ZScanSeedVertexFinderAlgorithm.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Vertexing/AdaptiveMultiVertexFinder.hpp"
#include "Acts/Vertexing/AdaptiveMultiVertexFitter.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/ImpactPointEstimator.hpp"
#include "Acts/Vertexing/LinearizedTrack.hpp"
#include "Acts/Vertexing/TrackDensityVertexFinder.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexFinderConcept.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "Acts/Vertexing/ZScanSeedVertexFinder.hpp"

#include "VertexingHelpers.hpp"

ActsExamples::ZScanSeedVertexFinderAlgorithm::ZScanSeedVertexFinderAlgorithm(
    const Config& cfg, Acts::Logging::Level lvl)
    : ActsExamples::BareAlgorithm("ZScanSeedVertexFinder", lvl), m_cfg(cfg) {
  if (m_cfg.inputSpacepoints.empty()) {
    throw std::invalid_argument(
        "You have to either provide seeds");
  }
  if (m_cfg.outputVertices.empty()) {
    throw std::invalid_argument("Missing output vertices collection");
  }
}

ActsExamples::ProcessCode ActsExamples::ZScanSeedVertexFinderAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // retrieve input seeds
  const auto& inputSpacepoints =
      ctx.eventStore.get<std::vector<SimSpacePoint>>(m_cfg.inputSpacepoints);

  // Setup the vertex fitter
  ZScanSeedVertexFinder<std::vector<SimSpacePoint>>::Config zscanSeedVtxCfg;
  ZScanSeedVertexFinder<std::vector<SimSpacePoint>> zscanSeedVertexFinder(zscanSeedVtxCfg);

  // find vertices and measure elapsed time
  auto t1 = std::chrono::high_resolution_clock::now();
  auto result = zscanSeedVertexFinder.findVertex(inputSpacepoints);
  auto t2 = std::chrono::high_resolution_clock::now();

  ACTS_INFO("Found " << result.size() << " vertices in event");
  for(auto r : result)
  {
    ACTS_INFO("Found vertex at z = " << r <<"mm");
  }

  // std::vector<Acts::Vertex<Acts::BoundTrackParameters>> vertices;
  // if (result.ok()) {
  //   vertices = std::move(result.value());
  // } else {
  //   ACTS_ERROR("Error in vertex finder: " << result.error().message());
  // }

  // show some debug output
  // ACTS_INFO("Found " << vertices.size() << " vertices in event");
  // for (const auto& vtx : vertices) {
  //   ACTS_INFO("Found vertex at " << vtx.fullPosition().transpose() << " with "
  //                                << vtx.tracks().size() << " tracks.");
  // }

  // store found vertices
  // ctx.eventStore.add(m_cfg.outputVertices, std::move(vertices));
  
  return ActsExamples::ProcessCode::SUCCESS;
}
