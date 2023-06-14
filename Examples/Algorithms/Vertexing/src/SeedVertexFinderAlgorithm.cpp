// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Vertexing/SeedVertexFinderAlgorithm.hpp"

#include "Acts/Vertexing/SeedVertexFinder.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <chrono>
#include <vector>

ActsExamples::SeedVertexFinderAlgorithm::SeedVertexFinderAlgorithm(
    const Config& cfg, Acts::Logging::Level lvl)
    : ActsExamples::IAlgorithm("SeedVertexFinder", lvl), m_cfg(cfg) {
  if (m_cfg.inputSpacepoints.empty()) {
    ACTS_ERROR("You have to either provide seeds");
  }
  if (m_cfg.outputVertices.empty()) {
    ACTS_ERROR("Missing output vertices collection");
  }

  m_inputSpacepoints.initialize(m_cfg.inputSpacepoints);
  m_outputSeedVertices.initialize(m_cfg.outputVertices);
}

ActsExamples::ProcessCode ActsExamples::SeedVertexFinderAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // retrieve input seeds
  const std::vector<ActsExamples::SimSpacePoint>& inputSpacepoints =
      m_inputSpacepoints(ctx);

  Acts::SeedVertexFinder<ActsExamples::SimSpacePoint>::Config zscanSeedVtxCfg;
  Acts::SeedVertexFinder<ActsExamples::SimSpacePoint> SeedVertexFinder(
      zscanSeedVtxCfg);

  // find vertices and measure elapsed time
  auto t1 = std::chrono::high_resolution_clock::now();
  auto vtx = SeedVertexFinder.findVertex(inputSpacepoints);
  auto t2 = std::chrono::high_resolution_clock::now();
  ACTS_INFO("Found a vertex in the event in " << (t2 - t1).count() / 1e6
                                              << " ms");
  ACTS_INFO("Found vertex at x = " << vtx[0] << "mm, y = " << vtx[1]
                                   << "mm, z = " << vtx[2] << "mm");

  std::vector<std::pair<Acts::Vector3, double>> results;
  results.push_back(std::make_pair(vtx, (t2 - t1).count() / 1e6));

  m_outputSeedVertices(ctx, std::move(results));

  return ActsExamples::ProcessCode::SUCCESS;
}
