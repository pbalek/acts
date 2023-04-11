// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/IAlgorithm.hpp"

#include <string>

namespace ActsExamples {

class ZScanSeedVertexFinderAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    /// Optional. Input spacepoints container.
    std::string inputSpacepoints;
    /// Output vertex collection
    std::string outputVertices;
  };

  ZScanSeedVertexFinderAlgorithm(const Config& cfg, Acts::Logging::Level lvl);

  /// @brief Find a vertex using spacepoints
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
