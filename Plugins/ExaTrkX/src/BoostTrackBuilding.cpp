// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/BoostTrackBuilding.hpp"

#include <map>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <torch/torch.h>

namespace {
template <typename vertex_t, typename edge_t, typename weight_t>
void weaklyConnectedComponents(vertex_t numNodes,
                               std::vector<vertex_t>& rowIndices,
                               std::vector<vertex_t>& colIndices,
                               std::vector<weight_t>& edgeWeights,
                               std::vector<vertex_t>& trackLabels) {
  typedef boost::adjacency_list<boost::vecS,         // edge list
                                boost::vecS,         // vertex list
                                boost::undirectedS,  // directedness
                                boost::no_property,  // property of vertices
                                float                // property of edges
                                >
      Graph;

  Graph g(numNodes);
  for (size_t idx = 0; idx < rowIndices.size(); ++idx) {
    boost::add_edge(rowIndices[idx], colIndices[idx], edgeWeights[idx], g);
  }

  [[maybe_unused]] size_t num_components =
      boost::connected_components(g, &trackLabels[0]);
}
}  // namespace

namespace Acts {

std::vector<std::vector<int>> BoostTrackBuilding::operator()(
    std::any, std::any edges, std::any weights,
    std::vector<int>& spacepointIDs) {
  ACTS_DEBUG("Start track building");
  const auto edgeTensor = std::any_cast<torch::Tensor>(edges).to(torch::kCPU);
  const auto edgeWeightTensor =
      std::any_cast<torch::Tensor>(weights).to(torch::kCPU);

  assert(edgeTensor.size(0) == 2);
  assert(edgeTensor.size(1) == edgeWeightTensor.size(0));

  const auto numSpacepoints = spacepointIDs.size();
  const auto numEdges = static_cast<std::size_t>(edgeWeightTensor.size(0));

  if (numEdges == 0) {
    ACTS_WARNING("No edges remained after edge classification");
    return {};
  }

  using vertex_t = int32_t;
  std::vector<vertex_t> rowIndices(numEdges);
  std::copy(edgeTensor.data_ptr<int64_t>(),
            edgeTensor.data_ptr<int64_t>() + numEdges, rowIndices.begin());

  std::vector<vertex_t> colIndices(numEdges);
  std::copy(edgeTensor.data_ptr<int64_t>() + numEdges,
            edgeTensor.data_ptr<int64_t>() + numEdges + numEdges,
            colIndices.begin());

  std::vector<float> edgeWeights(numEdges);
  std::copy(edgeWeightTensor.data_ptr<float>(),
            edgeWeightTensor.data_ptr<float>() + numEdges, edgeWeights.begin());

  std::vector<vertex_t> trackLabels(numSpacepoints);

  weaklyConnectedComponents<int32_t, int32_t, float>(
      numSpacepoints, rowIndices, colIndices, edgeWeights, trackLabels);

  ACTS_VERBOSE("Number of track labels: " << trackLabels.size());
  ACTS_VERBOSE("NUmber of unique track labels: " << [&]() {
    std::vector<vertex_t> sorted(trackLabels);
    std::sort(sorted.begin(), sorted.end());
    sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
    return sorted.size();
  }());

  if (trackLabels.size() == 0) {
    return {};
  }

  std::vector<std::vector<int>> trackCandidates;

  int existTrkIdx = 0;
  // map labeling from MCC to customized track id.
  std::map<int32_t, int32_t> trackLableToIds;

  for (auto idx = 0ul; idx < numSpacepoints; ++idx) {
    int32_t trackLabel = trackLabels[idx];
    int spacepointID = spacepointIDs[idx];

    int trkId;
    if (trackLableToIds.find(trackLabel) != trackLableToIds.end()) {
      trkId = trackLableToIds[trackLabel];
      trackCandidates[trkId].push_back(spacepointID);
    } else {
      // a new track, assign the track id
      // and create a vector
      trkId = existTrkIdx;
      trackCandidates.push_back(std::vector<int>{spacepointID});
      trackLableToIds[trackLabel] = trkId;
      existTrkIdx++;
    }
  }

  return trackCandidates;
}

}  // namespace Acts
