// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Vertexing/SeedVertexFinderAlgorithm.hpp"

#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "Acts/Vertexing/SeedVertexFinder.hpp"

#include <vector>
#include <chrono>

ActsExamples::SeedVertexFinderAlgorithm::SeedVertexFinderAlgorithm(
    const Config& cfg, Acts::Logging::Level lvl)
    : ActsExamples::IAlgorithm("SeedVertexFinder", lvl), m_cfg(cfg) {
  if (m_cfg.inputSpacepoints.empty()) { 
    ACTS_ERROR("You have to either provide seeds");
  }
  if (m_cfg.outputVertices.empty()) {
    ACTS_ERROR("Missing output vertices collection");
  }
}

ActsExamples::ProcessCode ActsExamples::SeedVertexFinderAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // retrieve input seeds
  const auto& inputSpacepoints = ctx.eventStore.get<SimSpacePointContainer>(m_cfg.inputSpacepoints);


  // Setup the seed vertex fitter, minimalization with respect to planes
  Acts::SeedVertexFinder<ActsExamples::SimSpacePoint>::Config zscanSeedVtxCfg1;
  Acts::SeedVertexFinder<ActsExamples::SimSpacePoint> SeedVertexFinder1(zscanSeedVtxCfg1);

  // find vertices and measure elapsed time
  auto t11 = std::chrono::high_resolution_clock::now();
  auto result1 = SeedVertexFinder1.findVertex(inputSpacepoints);
  auto t12 = std::chrono::high_resolution_clock::now();

  ACTS_INFO("Found a vertex in the event in "<<(t12-t11).count()/1e6<<" ms, using minimalization with respect to "<<zscanSeedVtxCfg1.minimalizeWRT);
  ACTS_INFO("Found vertex at x = " << result1[0] <<"mm, y = " << result1[1] <<"mm, z = " << result1[2] <<"mm");


  // Setup the seed vertex fitter, minimalization with respect to rays
  Acts::SeedVertexFinder<ActsExamples::SimSpacePoint>::Config zscanSeedVtxCfg2;
  zscanSeedVtxCfg2.minimalizeWRT="rays";
  Acts::SeedVertexFinder<ActsExamples::SimSpacePoint> SeedVertexFinder2(zscanSeedVtxCfg2);

  // find vertices and measure elapsed time
  auto t21 = std::chrono::high_resolution_clock::now();
  auto result2 = SeedVertexFinder2.findVertex(inputSpacepoints);
  auto t22 = std::chrono::high_resolution_clock::now();

  ACTS_INFO("Found a vertex in the event in "<<(t22-t21).count()/1e6<<" ms, using minimalization with respect to "<<zscanSeedVtxCfg2.minimalizeWRT);
  ACTS_INFO("Found vertex at x = " << result2[0] <<"mm, y = " << result2[1] <<"mm, z = " << result2[2] <<"mm");
  

  std::vector<Acts::Vector3> result{result1,result2};

  ctx.eventStore.add(m_cfg.outputVertices, std::move(result));
  
  return ActsExamples::ProcessCode::SUCCESS;
}
