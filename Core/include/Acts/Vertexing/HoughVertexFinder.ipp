// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cmath>
#include <system_error>
#include <fstream>

#include <Eigen/Eigenvalues>

template <typename spacepoint_t>
Acts::HoughVertexFinder<spacepoint_t>::HoughVertexFinder(
    const Acts::HoughVertexFinder<spacepoint_t>::Config& cfg,
    std::unique_ptr<const Logger> lgr)
    : m_cfg(cfg), m_logger(std::move(lgr)) {

  if(cfg.absEtaFractions.size() != cfg.absEtaRanges.size())
  {
    ACTS_ERROR("size of the absEtaFractions is "
              << cfg.absEtaFractions.size()
              << " but size of the absEtaRanges vector is "
              << cfg.absEtaRanges.size()
              <<"; these two have to be equal.");
  }

  if(cfg.rangeIterZ.size() != cfg.nBinsZIterZ.size())
  {
    ACTS_ERROR("size of the rangeIterZ is "
              << cfg.rangeIterZ.size()
              << " but size of the nBinsZIterZ vector is "
              << cfg.nBinsZIterZ.size()
              <<"; these two have to be equal.");
  }

  if(cfg.rangeIterZ.size() != cfg.nBinsCotThetaIterZ.size())
  {
    ACTS_ERROR("size of the rangeIterZ is "
              << cfg.rangeIterZ.size()
              << " but size of the nBinsCotThetaIterZ vector is "
              << cfg.nBinsCotThetaIterZ.size()
              <<"; these two have to be equal.");
  }
}

template <typename spacepoint_t>
Acts::Result<Acts::Vector3>
Acts::HoughVertexFinder<spacepoint_t>::find(
  const std::vector<spacepoint_t>& spacepoints) const {
  
  if(spacepoints.size()==0)
  {
    Acts::Result<Acts::Vector3>::failure(std::error_code());
  }

  Acts::ActsScalar absEtaRange=m_cfg.maxAbsEta;
  Acts::ActsScalar totalFrac=0.;
  for(unsigned int r=0;r<m_cfg.absEtaRanges.size();++r)
  {
    Acts::ActsScalar addToTotalFrac = m_cfg.absEtaFractions.at(r);
    if((totalFrac+addToTotalFrac)*spacepoints.size()>m_cfg.targetSPs)
    {
      Acts::ActsScalar needOnly = (m_cfg.targetSPs-totalFrac*spacepoints.size())/addToTotalFrac;
      absEtaRange = (r?m_cfg.absEtaRanges.at(r-1):0.) + (m_cfg.absEtaRanges.at(r)-(r?m_cfg.absEtaRanges.at(r-1):0.))*needOnly;
      break;
    }
    
    // Acts::ActsScalar thisRange = absEtaRanges.at(r) - (r?absEtaRanges.at(r-1):0.);
    // Acts::ActsScalar addToTotalFrac = absEtaFractions.at(r)*thisRange;
    // if((totalFrac+addToTotalFrac)*spacepoints.size()>m_cfg.targetSPs)
    // {
    //   Acts::ActsScalar needOnly = (m_cfg.targetSPs-totalFrac*spacepoints.size())/addToTotalFrac; // bug!!!
    //   absEtaRange = (r?absEtaRanges.at(r-1):0.) + needOnly*thisRange;
    //   break;
    // }

    totalFrac += addToTotalFrac;
  }
  if(absEtaRange>m_cfg.maxAbsEta) absEtaRange = m_cfg.maxAbsEta;
  if(absEtaRange<m_cfg.minAbsEta) absEtaRange = m_cfg.minAbsEta;

  const Acts::ActsScalar maxCotTheta = std::sinh(absEtaRange);
  const Acts::ActsScalar minCotTheta = -1.*maxCotTheta;

  Acts::ActsScalar binsNumDecrease=1.;
  if(spacepoints.size()*totalFrac < m_cfg.targetSPs)
  {
    binsNumDecrease = std::pow(m_cfg.binsCotThetaDecrease, std::log(m_cfg.targetSPs/(spacepoints.size()*totalFrac)));
  }

  Acts::Vector3 vtx = m_cfg.defVtxPosition;

  for(std::uint32_t iter=0;iter<m_cfg.rangeIterZ.size();++iter)
  {
    auto vtxNew = findHoughPeak(spacepoints,
                                vtx, m_cfg.rangeIterZ.at(iter), m_cfg.nBinsZIterZ.at(iter),
                                minCotTheta, maxCotTheta, m_cfg.nBinsCotThetaIterZ.at(iter)/binsNumDecrease);
    
    if(!vtxNew.ok())
    {
      // vertex not found
      Acts::Result<Acts::Vector3>::failure(std::error_code());
    }

    vtx = vtxNew.value();
  }

  return Acts::Result<Acts::Vector3>::success(vtx);
}


template <typename spacepoint_t>
Acts::Result<Acts::Vector3> Acts::HoughVertexFinder<spacepoint_t>::findHoughPeak(
    const std::vector<spacepoint_t>& spacepoints,
    Acts::Vector3 vtxOld, Acts::ActsScalar rangeZ, std::uint32_t numZBins,
    Acts::ActsScalar minCotTheta, Acts::ActsScalar maxCotTheta, std::uint32_t numCotThetaBins) const {

  const Acts::ActsScalar zBinSize = 2.*rangeZ/numZBins;
  const Acts::ActsScalar invCotThetaBinSize = numCotThetaBins/(maxCotTheta-minCotTheta);
  const Acts::ActsScalar minZ=vtxOld[2]-rangeZ, maxZ=vtxOld[2]+rangeZ;

  MultiIndexedVector2D<std::uint32_t> houghImage(numZBins,numCotThetaBins,0);

  std::vector<std::uint32_t> houghZProjection(numZBins,0);

  std::vector<Acts::ActsScalar> vtxZPositions;
  for(std::uint32_t zBin=0;zBin<numZBins;zBin++) {
    vtxZPositions.push_back(Acts::HoughTransformUtils::binCenter(
      minZ, maxZ, numZBins, zBin
    ) );
  }

  for (const auto& sp : spacepoints) {
    if(sp.z()>maxZ)
    {
      if((sp.z()-maxZ)/sp.r() > maxCotTheta) {
        continue;
      }
    }
    else if(sp.z()<minZ)
    {
      if((sp.z()-minZ)/sp.r() < minCotTheta) {
        continue;
      }
    }

    Acts::ActsScalar sp_invr=1./std::sqrt((sp.x()-vtxOld[0])*(sp.x()-vtxOld[0]) + (sp.y()-vtxOld[1])*(sp.y()-vtxOld[1]));

    std::uint32_t zFrom = static_cast<std::uint32_t>(((sp.z()-maxCotTheta/sp_invr)-minZ)/zBinSize +1);  
    // if(zFrom<0) zFrom=0; // TODO: it really can't happen??
    std::uint32_t zTo   = static_cast<std::uint32_t>(((sp.z()-minCotTheta/sp_invr)-minZ)/zBinSize);
    if(zTo>=numZBins) zTo=numZBins;

    for(std::uint32_t zBin=zFrom;zBin<zTo;zBin++) {
      Acts::ActsScalar cotTheta = (sp.z()-vtxZPositions[zBin])*sp_invr;

      std::uint32_t cotThetaBin = static_cast<std::uint32_t>((cotTheta-minCotTheta)*invCotThetaBinSize);

      std::uint32_t cotThetaFrom = std::max<std::uint32_t>(cotThetaBin-m_cfg.fillNeighbours,0);
      std::uint32_t cotThetaTo   = std::min(cotThetaBin+m_cfg.fillNeighbours+1,numCotThetaBins);

      for(std::uint32_t cotBin=cotThetaFrom; cotBin<cotThetaTo; ++cotBin)
      {
        ++houghImage(zBin,cotBin);
      }
    }
  }

  std::uint32_t maxZBin=0;
  for(std::uint32_t zBin=0;zBin<numZBins;zBin++) {
    houghZProjection[zBin] = std::reduce(houghImage[zBin], houghImage[zBin]+numCotThetaBins, static_cast<std::uint32_t>(0),
    [&](std::uint32_t& lhs, std::uint32_t& rhs) { return (lhs>=m_cfg.minHits)*lhs + (rhs>=m_cfg.minHits)*rhs;} );

    if(houghZProjection[zBin] > houghZProjection.at(maxZBin)) maxZBin=zBin;
  }

  Acts::ActsScalar avg=std::accumulate(houghZProjection.begin(), houghZProjection.end(), 0.)/houghZProjection.size();
  Acts::ActsScalar sumEntries=0;
  Acts::ActsScalar meanZPeak=0.;

  for(std::uint32_t zBin=std::max<std::uint32_t>(maxZBin-m_cfg.peakWidth,0); 
      zBin<=std::min(numZBins-1,maxZBin+m_cfg.peakWidth); 
      ++zBin) {
    sumEntries += std::max(houghZProjection.at(zBin)-avg, 0.);
    meanZPeak   += vtxZPositions[zBin]*std::max<std::uint32_t>(houghZProjection.at(zBin)-avg, 0.);
  }

  if(sumEntries!=0.)
  {
    meanZPeak/=sumEntries;
    return Acts::Result<Acts::Vector3>::success({vtxOld[0], vtxOld[1], meanZPeak});
  }

  // vertex not found; Hough image empty
  return Acts::Result<Acts::Vector3>::failure(std::error_code());
}
