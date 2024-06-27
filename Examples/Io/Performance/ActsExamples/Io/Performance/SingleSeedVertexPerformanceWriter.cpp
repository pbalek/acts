// This file is part of the Acts project.
//
// Copyright (C) 2019-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Performance/SingleSeedVertexPerformanceWriter.hpp"

// #include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/AverageSimHits.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"

#include <ios>
#include <map>
#include <set>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>

ActsExamples::SingleSeedVertexPerformanceWriter::SingleSeedVertexPerformanceWriter(
    const ActsExamples::SingleSeedVertexPerformanceWriter::Config& config,
    Acts::Logging::Level level)
    : WriterT(config.inputVertices, "SingleSeedVertexPerformanceWriter", level),
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
    throw std::invalid_argument(
        "Collection with selected truth particles missing");
  }
  if (m_cfg.inputVertices.empty()) {
    throw std::invalid_argument("Collection with vertices missing");
  }

  m_inputAllTruthParticles.initialize(m_cfg.inputAllTruthParticles);
  m_inputSelectedTruthParticles.initialize(m_cfg.inputSelectedTruthParticles);
  if (!m_cfg.inputTime.empty()) {
    m_inputTime.initialize(m_cfg.inputTime);
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
    m_outputTree->Branch("diffx", &m_diffX);
    m_outputTree->Branch("diffy", &m_diffY);
    m_outputTree->Branch("diffz", &m_diffZ);

    m_outputTree->Branch("recoX", &m_recoX);
    m_outputTree->Branch("recoY", &m_recoY);
    m_outputTree->Branch("recoZ", &m_recoZ);

    m_outputTree->Branch("truthX", &m_truthX);
    m_outputTree->Branch("truthY", &m_truthY);
    m_outputTree->Branch("truthZ", &m_truthZ);

    // for(int i = 0; i <100;++i)
    // {
    //   m_outputTree->Branch(Form("distance_i%d",i), &m_distance[i]);
    // }
    // m_outputTree->Branch("distanceSize", &m_distanceSize);
    // m_outputTree->Branch("diffx_i", &m_diffX_i);
    // m_outputTree->Branch("diffy_i", &m_diffY_i);
    // m_outputTree->Branch("diffz_i", &m_diffZ_i);

    m_outputTree->Branch("nRecoVtx", &m_nrecoVtx);
    m_outputTree->Branch("nTrueVtx", &m_ntrueVtx);
    m_outputTree->Branch("nVtxDetectorAcceptance", &m_nVtxDetAcceptance);
    m_outputTree->Branch("timeMS", &m_timeMS);

    m_outputTree->Branch("SPNum", &m_SPNum);
    m_outputTree->Branch("usedSPNum", &m_usedSPNum);
    // m_outputTree->Branch("farSPNum", &m_farSPNum);
    // m_outputTree->Branch("tripletNum", &m_tripletNum);
  }
}

ActsExamples::SingleSeedVertexPerformanceWriter::
    ~SingleSeedVertexPerformanceWriter() {
  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
}

ActsExamples::ProcessCode
ActsExamples::SingleSeedVertexPerformanceWriter::finalize() {
  m_outputFile->cd();
  m_outputTree->Write();
  m_outputFile->Close();

  return ProcessCode::SUCCESS;
}

int ActsExamples::SingleSeedVertexPerformanceWriter::getNumberOfTruePriVertices(
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

ActsExamples::ProcessCode ActsExamples::SingleSeedVertexPerformanceWriter::writeT(
    const AlgorithmContext& ctx,
    const std::vector<std::vector<double>>& vertices) {
  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  m_nrecoVtx = vertices.size();
  ACTS_DEBUG("Number of reco vertices in event: " << m_nrecoVtx);

  // // Read truth particle input collection
  const auto& allTruthParticles = m_inputAllTruthParticles(ctx);
  // Get number of generated true primary vertices
  m_ntrueVtx = getNumberOfTruePriVertices(allTruthParticles);
  ACTS_VERBOSE("Total number of generated truth particles in event : "
               << allTruthParticles.size());

  // Read selected truth particle input collection
  const auto& selectedTruthParticles = m_inputSelectedTruthParticles(ctx);
  // // Get number of detector-accepted true primary vertices
  m_nVtxDetAcceptance = getNumberOfTruePriVertices(selectedTruthParticles);
  ACTS_VERBOSE("Total number of detector-accepted truth primary vertices : "
               << m_nVtxDetAcceptance);

  // Loop over truth particles and find the truth
  std::vector<int> truthVertices;
  for (const auto& particle : allTruthParticles) {
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

  // std::cout<<"vector vertices.size = "<<vertices.size()<<std::endl;
  std::vector<std::vector<double>> verticesReal={vertices.at(0)};

  std::vector<double> time;
  time.push_back(vertices.at(1).at(0));
  // std::cout<<"time, size = "<< time.size()<<" : ";
  // for(unsigned int j=0;j<time.size();++j)
  // {
  //   std::cout<<time.at(j)<<" ";
  // }
  // std::cout<<std::endl;

  // std::vector<double> numSPs=vertices.at(102);
  // std::cout<<"numSPs, size = "<< numSPs.size()<<" : ";
  // for(unsigned int j=0;j<numSPs.size();++j)
  // {
  //   std::cout<<numSPs.at(j)<<" ";
  // }
  // std::cout<<std::endl;

  // for(unsigned int i=0;i<vertices.size();++i)
  // {
  //   std::cout<<"vertex "<<i<<", vector size "<<vertices.at(i).size()<<std::endl;
  //   std::cout<<"    ";
  //   for(unsigned int j=0;j<vertices.at(i).size();++j)
  //   {
  //     std::cout<<vertices.at(i).at(j)<<" ";
  //   }
  //   std::cout<<std::endl;
  //   if(i>=200) break;
  // }

  for (const auto& vtx : verticesReal) {
    // find a particle with such truth vertex

    for (const auto& particle : allTruthParticles) {
      int priVtxId = particle.particleId().vertexPrimary();

      if (priVtxId != maxOccurrenceId) {
        continue;
      }
      // this is the truth vertex
      const auto& truePos = particle.position();

      m_diffX.push_back(vtx[0] - truePos[0]);
      m_diffY.push_back(vtx[1] - truePos[1]);
      m_diffZ.push_back(vtx[2] - truePos[2]);

      m_truthX.push_back(truePos[0]);
      m_truthY.push_back(truePos[1]);
      m_truthZ.push_back(truePos[2]);

      m_recoX.push_back(vtx[0]);
      m_recoY.push_back(vtx[1]);
      m_recoZ.push_back(vtx[2]);

      m_timeMS.push_back(time[0]);

      m_SPNum.push_back(vertices.at(1).at(1));
      m_usedSPNum.push_back(vertices.at(1).at(2));
      // m_farSPNum.push_back(numSPs[2]);
      // m_tripletNum.push_back(numSPs[3]);


      // for(unsigned int tr=2;tr<102;tr++)
      // {
      //   // std::cout<<" diff i, tr "<<tr<<std::endl;
      //   auto vtx_i=vertices.at(tr);
      //   m_diffX_i.push_back(vtx_i[0] - truePos[0]);
      //   m_diffY_i.push_back(vtx_i[1] - truePos[1]);
      //   m_diffZ_i.push_back(vtx_i[2] - truePos[2]);
      // }

      // for(unsigned int tr=103;tr<vertices.size();tr++)
      // {
      //   // std::cout<<" triplet, tr "<<tr<<std::endl;
      //   auto triplet=vertices.at(tr);

      //   int iter=(int)(triplet[0]+0.01);
      //   double * abg = &(triplet[1]);
      //   double * delta = &(triplet[4]);
      //   double * startPoint = &(triplet[5]);
      //   double * direction  = &(triplet[8]);

      //   double effectEccSq=(triplet[0]-iter)*10.0;

      //   double abg_dot_vtx = abg[0]*truePos[0]+abg[1]*truePos[1]+abg[2]*truePos[2] + delta[0];

      //   double vtxDiff[3]={truePos[0]-startPoint[0], truePos[1]-startPoint[1], truePos[2]-startPoint[2]};
      //   double diff_cross[3]={vtxDiff[1]*direction[2]-vtxDiff[2]*direction[1],
      //                         vtxDiff[2]*direction[0]-vtxDiff[0]*direction[2],
      //                         vtxDiff[0]*direction[1]-vtxDiff[1]*direction[0]};

      //   double distanceSq=effectEccSq*abg_dot_vtx*abg_dot_vtx + (1.-effectEccSq)*(diff_cross[0]*diff_cross[0]+diff_cross[1]*diff_cross[1]+diff_cross[2]*diff_cross[2]);

      //   // ACTS_INFO(iter<<" "<<tr<<" : abg = "<<abg[0]<<", "<<abg[1]<<", "<<abg[2]<<"; delta = "<<delta<<"; startPoint "<<startPoint[0]<<", "<<startPoint[1]<<", "<<startPoint[2]<<", direction = "<<direction[0]<<", "<<direction[1]<<", "<<direction[2]);
      //   // ACTS_INFO(iter<<" "<<tr<<" : effectEccSq = "<<effectEccSq<<", dist sq to plane "<<abg_dot_vtx * abg_dot_vtx <<", dist sq to ray "<< (diff_cross[0]*diff_cross[0]+diff_cross[1]*diff_cross[1]+diff_cross[2]*diff_cross[2]));

      //   m_distance[iter].push_back(std::sqrt(distanceSq));
      // }

      // for(unsigned int iter=0;iter<100;iter++)
      // {
      //   m_distanceSize.push_back(m_distance[iter].size());
      // }

      break;  // particle loop
    }
  }

  // Retrieve and set reconstruction time
  // if (!m_cfg.inputTime.empty()) {
  //   const auto& reconstructionTimeMS = m_inputTime(ctx);
  //   m_timeMS = reconstructionTimeMS;
  // } else {
  //   m_timeMS = -1;
  // }

  // fill the variables
  m_outputTree->Fill();

  m_diffX.clear();
  m_diffY.clear();
  m_diffZ.clear();
  m_truthX.clear();
  m_truthY.clear();
  m_truthZ.clear();
  m_recoX.clear();
  m_recoY.clear();
  m_recoZ.clear();
  // for(int i = 0; i <100;++i)
  // {
  //   m_distance[i].clear();
  // }
  // m_distanceSize.clear();
  // m_diffX_i.clear();
  // m_diffY_i.clear();
  // m_diffZ_i.clear();
  m_timeMS.clear();
  m_SPNum.clear();
  m_usedSPNum.clear();
  // m_farSPNum.clear();
  // m_tripletNum.clear();

  // std::cout<<"writing done."<<std::endl;

  return ProcessCode::SUCCESS;
} 