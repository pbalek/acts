// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Vertexing/Vertex.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <mutex>
#include <vector>

class TFile;
class TTree;

namespace ActsExamples {

/// @class SingleSeedVertexPerformanceWriter
///
/// Writes out the number of reconstructed seed vertices along with
/// the number of truth vertices in detector acceptance.
/// Additionally it matches the reco vertices to their truth vertices
/// and write out the difference in x, y, and z position.
class SingleSeedVertexPerformanceWriter final
    : public WriterT<std::vector<std::vector<double>>> {
 public:
  struct Config {
    /// All input truth particle collection.
    std::string inputAllTruthParticles;
    /// Selected input truth particle collection.
    std::string inputSelectedTruthParticles;
    /// Input vertex collection.
    std::string inputVertices;
    /// Input reconstruction time.
    std::string inputTime;
    /// Output filename.
    std::string filePath = "seedvertexingperformance.root";
    /// Name of the output tree.
    std::string treeName = "seedvertextree";
    /// File access mode.
    std::string fileMode = "RECREATE";
  };

  /// Constructor
  ///
  /// @param config Configuration struct
  /// @param level Message level declaration
  SingleSeedVertexPerformanceWriter(const Config& config,
                                  Acts::Logging::Level level);

  ~SingleSeedVertexPerformanceWriter() override;

  /// End-of-run hook
  ProcessCode finalize() override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 protected:
  /// @brief Write method called by the base class
  /// @param [in] ctx is the algorithm context for event information
  ProcessCode writeT(
      const AlgorithmContext& ctx,
      const std::vector<std::vector<double>>& vertices) override;

 private:
  Config m_cfg;             ///< The config class
  std::mutex m_writeMutex;  ///< Mutex used to protect multi-threaded writes
  TFile* m_outputFile{nullptr};  ///< The output file
  TTree* m_outputTree{nullptr};  ///< The output tree

  std::vector<float>
      m_diffX;  ///< Difference in x position between reco and true vtx
  std::vector<float>
      m_diffY;  ///< Difference in y position between reco and true vtx
  std::vector<float>
      m_diffZ;  ///< Difference in z position between reco and true vtx

  std::vector<float> m_truthX;
  std::vector<float> m_truthY;
  std::vector<float> m_truthZ;

  std::vector<float> m_recoX;
  std::vector<float> m_recoY;
  std::vector<float> m_recoZ;

  int m_nrecoVtx = -1;           ///< Number of reconstructed vertices
  int m_ntrueVtx = -1;           ///< Number of true vertices
  int m_nVtxDetAcceptance = -1;  ///< Number of vertices in detector acceptance
  std::vector<double> m_timeMS;  ///< Reconstruction time in ms
  std::vector<double> m_SPNum;
  std::vector<double> m_usedSPNum;
  std::vector<double> m_peakSPNum;
  // std::vector<double> m_tripletNum;
  // std::vector<double> m_distance[100];
  // std::vector<double> m_distanceSize;

  // std::vector<float> m_diffX_i, m_diffY_i, m_diffZ_i;

  int getNumberOfTruePriVertices(const SimParticleContainer& collection) const;

  ReadDataHandle<SimParticleContainer> m_inputAllTruthParticles{
      this, "particles_input"};
  ReadDataHandle<SimParticleContainer> m_inputSelectedTruthParticles{
      this, "particles_selected"};
  ReadDataHandle<int> m_inputTime{this, ""};
};

}  // namespace ActsExamples