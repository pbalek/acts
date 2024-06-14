// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Vertexing/SingleSeedVertexFinderAlgorithm.hpp"

#include "Acts/Vertexing/SingleSeedVertexFinder.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <chrono>
#include <vector>

ActsExamples::SingleSeedVertexFinderAlgorithm::SingleSeedVertexFinderAlgorithm(
    const Config& cfg, Acts::Logging::Level lvl)
    : ActsExamples::IAlgorithm("SingleSeedVertexFinder", lvl), m_cfg(cfg) {
  if (m_cfg.inputSpacepoints.empty()) {
    ACTS_ERROR("You have to provide seeds");
  }
  if (m_cfg.outputVertices.empty()) {
    ACTS_ERROR("Missing output vertices collection");
  }

  m_inputSpacepoints.initialize(m_cfg.inputSpacepoints);
  // m_outputVertices.initialize(m_cfg.outputVertices);
  m_outputSeedVertices.initialize(m_cfg.outputVertices);
}

ActsExamples::ProcessCode
ActsExamples::SingleSeedVertexFinderAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  

  // FILE* file = fopen("/proc/self/status", "r");
  // char line[128];
  // while (fgets(line, 128, file) != NULL){
  //     if (strncmp(line, "VmSize:", 7) == 0){
  //         ACTS_INFO("Beginning of SingleSeedVertexFinderAlgorithm::execute: VmSize="<<line);
  //     }
  //     if (strncmp(line, "VmRSS:", 6) == 0){
  //       ACTS_INFO("Beginning of SingleSeedVertexFinderAlgorithm::execute: VmRSS="<<line);
  //     }
  // }
  // fclose(file);

  // retrieve input seeds
  const std::vector<ActsExamples::SimSpacePoint>& inputSpacepoints =
      m_inputSpacepoints(ctx);

  Acts::SingleSeedVertexFinder<ActsExamples::SimSpacePoint>::Config
      singleSeedVtxCfg;
  singleSeedVtxCfg.mixedEccentricity=std::sqrt(1.-1./(m_cfg.ecc*m_cfg.ecc));

  Acts::SingleSeedVertexFinder<ActsExamples::SimSpacePoint>
      SingleSeedVertexFinder(singleSeedVtxCfg);

  // find vertices and measure elapsed time
  // file = fopen("/proc/self/status", "r");
  // while (fgets(line, 128, file) != NULL){
  //     if (strncmp(line, "VmSize:", 7) == 0){
  //         ACTS_INFO("Before SingleSeedVertexFinder: VmSize="<<line);
  //     }
  //     if (strncmp(line, "VmRSS:", 6) == 0){
  //       ACTS_INFO("Before SingleSeedVertexFinder: VmRSS="<<line);
  //     }
  // }
  // fclose(file);

  auto t1 = std::chrono::high_resolution_clock::now();
  auto [vtx, rejectVector] = SingleSeedVertexFinder.findVertex(inputSpacepoints);
  auto t2 = std::chrono::high_resolution_clock::now();
  // file = fopen("/proc/self/status", "r");
  // while (fgets(line, 128, file) != NULL){
  //     if (strncmp(line, "VmSize:", 7) == 0){
  //         ACTS_INFO("After SingleSeedVertexFinder: VmSize="<<line);
  //     }
  //     if (strncmp(line, "VmRSS:", 6) == 0){
  //       ACTS_INFO("After SingleSeedVertexFinder: VmRSS="<<line);
  //     }
  // }
  // fclose(file);

  ACTS_INFO("rejectVector "<<rejectVector.size());

  if (vtx.ok()) {
    ACTS_INFO("Found a vertex in the event in " << (t2 - t1).count() / 1e6
                                                << " ms");
    ACTS_INFO("Found vertex at x = " << vtx.value()[0]
                                     << "mm, y = " << vtx.value()[1]
                                     << "mm, z = " << vtx.value()[2] << "mm");

    // std::vector<Acts::Vertex<Acts::BoundTrackParameters>> vertexCollection;
    // vertexCollection.emplace_back(vtx.value());

    // std::vector<std::pair<Acts::Vector3, double>> results;
    // results.push_back(std::make_pair(vtx.value(), (t2 - t1).count() / 1e6));

    std::vector<std::vector<double>> results;
    results.push_back({vtx.value()[0], vtx.value()[1], vtx.value()[2]});
    results.push_back({(t2 - t1).count() / 1e6});
    for(auto trip : rejectVector)
    {
      results.push_back(trip);
    }

    m_outputSeedVertices(ctx, std::move(results));

    // store found vertices
    // m_outputVertices(ctx, std::move(vertexCollection));
  } else {
    ACTS_INFO("Not found a vertex in the event after "
              << (t2 - t1).count() / 1e6 << " ms");
    
    // auto vtx_dummy = Acts::Result<Acts::Vector3>::success(Acts::Vector3::Zero());
    // std::vector<Acts::Vertex<Acts::BoundTrackParameters>> vertexCollection;
    // vertexCollection.emplace_back(vtx_dummy.value());
    // m_outputVertices(ctx, std::move(vertexCollection));

    // std::vector<std::pair<Acts::Vector3, double>> results{{{0.,0.,0.},-1}};
    std::vector<std::vector<double>> results{{0.,0.,0.},{-1}};
    for(auto trip : rejectVector)
    {
      results.push_back(trip);
    }

    m_outputSeedVertices(ctx, std::move(results));
  }

  return ActsExamples::ProcessCode::SUCCESS;
}
