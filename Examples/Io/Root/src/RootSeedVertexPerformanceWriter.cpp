// This file is part of the Acts project.
//
// Copyright (C) 2019-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootSeedVertexPerformanceWriter.hpp"

// #include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/AverageSimHits.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"

#include <ios>
#include <stdexcept>
#include <set>
#include <map>

#include <TFile.h>
#include <TTree.h>

ActsExamples::RootSeedVertexPerformanceWriter::RootSeedVertexPerformanceWriter(
    const ActsExamples::RootSeedVertexPerformanceWriter::Config& config,
    Acts::Logging::Level level)
    : WriterT(config.inputVertices, "RootSeedVertexPerformanceWriter", level),
      m_cfg(config) {
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing output filename");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }
  if (m_cfg.inputAllTruthParticles.empty()) {
    throw std::invalid_argument("Collection with all truth particles missing");
  }
  if (m_cfg.inputSelectedTruthParticles.empty()) {
    throw std::invalid_argument("Collection with selected truth particles missing");
  }
  if (m_cfg.inputVertices.empty()) {
    throw std::invalid_argument("Collection with vertices missing");
  }

  // Setup ROOT I/O
  auto path = m_cfg.filePath;
  m_outputFile = TFile::Open(path.c_str(), m_cfg.fileMode.c_str());
  if (m_outputFile == nullptr) {
    throw std::ios_base::failure("Could not open '" + path);
  }
  m_outputFile->cd();
  m_outputTree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());
  if (m_outputTree == nullptr) {
    throw std::bad_alloc();
  } else {
    // I/O parameters
    m_outputTree->Branch("diffx", &m_diffx);
    m_outputTree->Branch("diffy", &m_diffy);
    m_outputTree->Branch("diffz", &m_diffz);

    m_outputTree->Branch("recoX", &m_recoX);
    m_outputTree->Branch("recoY", &m_recoY);
    m_outputTree->Branch("recoZ", &m_recoZ);

    m_outputTree->Branch("truthX", &m_truthX);
    m_outputTree->Branch("truthY", &m_truthY);
    m_outputTree->Branch("truthZ", &m_truthZ);

    m_outputTree->Branch("nRecoVtx", &m_nrecoVtx);
    m_outputTree->Branch("nTrueVtx", &m_ntrueVtx);
    m_outputTree->Branch("nVtxDetectorAcceptance", &m_nVtxDetAcceptance);
    m_outputTree->Branch("timeMS", &m_timeMS);
  }
}

ActsExamples::RootSeedVertexPerformanceWriter::~RootSeedVertexPerformanceWriter() {
  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
}

ActsExamples::ProcessCode
ActsExamples::RootSeedVertexPerformanceWriter::finalize() {
  m_outputFile->cd();
  m_outputTree->Write();
  m_outputFile->Close();

  return ProcessCode::SUCCESS;
}


int ActsExamples::RootSeedVertexPerformanceWriter::getNumberOfTruePriVertices(
    const SimParticleContainer& collection) const {
  // Vector to store indices of all primary vertices
  std::set<int> allPriVtxIds;
  for (const auto& p : collection) {
    int priVtxId = p.particleId().vertexPrimary();
    int secVtxId = p.particleId().vertexSecondary();
    if (secVtxId != 0) {
      // truthparticle from secondary vtx
      continue;
    }
    // Insert to set, removing duplicates
    allPriVtxIds.insert(priVtxId);
  }
  // Size of set corresponds to total number of primary vertices
  return allPriVtxIds.size();
}

ActsExamples::ProcessCode ActsExamples::RootSeedVertexPerformanceWriter::writeT(
    const AlgorithmContext& ctx,
    const std::vector<Acts::Vector3>& vertices) {
  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  m_nrecoVtx = vertices.size();
  ACTS_DEBUG("Number of reco vertices in event: " << m_nrecoVtx);

  // // Read truth particle input collection
  const auto& allTruthParticles = ctx.eventStore.get<SimParticleContainer>(m_cfg.inputAllTruthParticles);
  // Get number of generated true primary vertices
  m_ntrueVtx = getNumberOfTruePriVertices(allTruthParticles);
  ACTS_VERBOSE("Total number of generated truth particles in event : " << allTruthParticles.size());

  // Read selected truth particle input collection
  const auto& selectedTruthParticles = ctx.eventStore.get<SimParticleContainer>(
      m_cfg.inputSelectedTruthParticles);
  // // Get number of detector-accepted true primary vertices
  m_nVtxDetAcceptance = getNumberOfTruePriVertices(selectedTruthParticles);
  ACTS_VERBOSE("Total number of detector-accepted truth primary vertices : " << m_nVtxDetAcceptance);

  // Loop over truth particles and find the truth
  std::vector<int> truthVertices;
  for (const auto& particle : allTruthParticles)
  {
    int priVtxId = particle.particleId().vertexPrimary();
    truthVertices.push_back(priVtxId);
  }
  // Find truth vertex with most particles
  std::map<int, int> fmap;
  for (int priVtxId : truthVertices) {
    fmap[priVtxId]++;
  }
  int maxOccurrence = -1;
  int maxOccurrenceId = -1;
  for (auto it : fmap) {
    if (it.second > maxOccurrence) {
      maxOccurrence = it.second;
      maxOccurrenceId = it.first;
    }
  }

  // Loop over all reco vertices and find associated truth particles
  std::vector<SimParticleContainer> truthParticlesAtVtxContainer;
  for (const auto& vtx : vertices) 
  {
    // find a particle with such truth vertex
    for (const auto& particle : allTruthParticles)
    {
      int priVtxId = particle.particleId().vertexPrimary();

      if (priVtxId != maxOccurrenceId) continue;
      // this is the truth vertex
      const auto& truePos = particle.position();

      m_diffx.push_back(vtx[0] - truePos[0]);
      m_diffy.push_back(vtx[1] - truePos[1]);
      m_diffz.push_back(vtx[2] - truePos[2]);

      m_truthX.push_back(truePos[0]);
      m_truthY.push_back(truePos[1]);
      m_truthZ.push_back(truePos[2]);

      m_recoX.push_back(vtx[0]);
      m_recoY.push_back(vtx[1]);
      m_recoZ.push_back(vtx[2]);

      break; // particle loop
    }
  }

  // Retrieve and set reconstruction time
  if (!m_cfg.inputTime.empty())
  {
    const auto& reconstructionTimeMS = ctx.eventStore.get<int>(m_cfg.inputTime);
    m_timeMS = reconstructionTimeMS;
  }
  else
  {
    m_timeMS = -1;
  }

  // fill the variables
  m_outputTree->Fill();

  m_diffx.clear();
  m_diffy.clear();
  m_diffz.clear();
  m_truthX.clear();
  m_truthY.clear();
  m_truthZ.clear();
  m_recoX.clear();
  m_recoY.clear();
  m_recoZ.clear();

  return ProcessCode::SUCCESS;
}
