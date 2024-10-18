// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cmath>
#include <system_error>
#include <fstream>

#include <Eigen/Eigenvalues>

template <typename spacepoint_t>
Acts::SingleSeedVertexFinder<spacepoint_t>::SingleSeedVertexFinder(
    const Acts::SingleSeedVertexFinder<spacepoint_t>::Config& cfg,
    std::unique_ptr<const Logger> lgr)
    : m_cfg(cfg), m_logger(std::move(lgr)) {
  if (m_cfg.minNumTriplets<100) {
    ACTS_INFO("value of minNumTriplets is "
              << cfg.minNumTriplets
              << "; this algorithm needs at least few hundreds triplets to work properly.");
  }
  if (cfg.numPhiSlices < 3) {
    ACTS_INFO("value of numPhiSlices is "
              << cfg.numPhiSlices
              << ", which is less than 3. There will be duplicate triplets.");
  }
  if (cfg.useFracPhiSlices <= 0. || cfg.useFracPhiSlices > 1.) {
    ACTS_ERROR("value of useFracPhiSlices is "
               << cfg.useFracPhiSlices
               << ", allowed values are between 0 and 1");
  }
  if (cfg.useFracZSlices <= 0. || cfg.useFracZSlices > 1.) {
    ACTS_ERROR("value of useFracZSlices is "
               << cfg.useFracZSlices << ", allowed values are between 0 and 1");
  }
  if(cfg.minimalizeWRT != "planes" && cfg.minimalizeWRT != "rays" && cfg.minimalizeWRT != "mixed" && cfg.minimalizeWRT != "hough") {
    ACTS_ERROR("value of minimalizeWRT is "
               << cfg.minimalizeWRT
               << ", allowed values are \"planes\" or \"rays\" or \"mixed\" or \"hough\".");
  }
  if(cfg.minimalizeWRT == "mixed")
  {
    if(cfg.mixedEccentricity < 0.0 || cfg.mixedEccentricity > 1.0) {
      ACTS_ERROR("value of mixedEccentricity is "
              << cfg.mixedEccentricity
              << ", allowed values are between 0 and 1.");
    }

    m_effectEccSq = cfg.mixedEccentricity*cfg.mixedEccentricity;
  }
  else {
    if(cfg.minimalizeWRT != "hough")
    {
      m_effectEccSq = (cfg.minimalizeWRT=="planes" ? 1. : 0.);
    }
  }
  if ((cfg.removeFraction < 0. || cfg.removeFraction >= 1.) && cfg.minimalizeWRT != "hough") {
    ACTS_ERROR("value of removeFraction is "
               << cfg.removeFraction << ", allowed values are between 0 and 1");
  }
}

template <typename spacepoint_t>
// std::pair<Acts::Result<Acts::Vector3>, std::vector<std::vector<Acts::ActsScalar>>>
Acts::Result<std::vector<double>>
Acts::SingleSeedVertexFinder<spacepoint_t>::findVertex(
    const std::vector<spacepoint_t>& spacepoints) const {
  // ACTS_INFO("Have "<<spacepoints.size()<<" spacepoints, size = "<<spacepoints.size()*sizeof(spacepoints[0])<<"B");

  // minTheta=1.0: eta=0.6;  minTheta=0.7: eta=1.0;  minTheta=0.27: eta=2.0; 
  //   usedSP=SP*0.162    ;    usedSP=SP*0.217    ;    usedSP=SP*0.435;
  // cca. 0.2*SP per unit of |eta| will be usedSP, raise to 0.3*SP for |eta|<2.0

  // if(spacepoints.size()<2000 && m_cfg.mixedEccentricity<0.5)
  // {
  //   char outname[50];
  //   sprintf(outname, "hough_transform_SPs_%ld.txt", spacepoints.size());
  //   std::ofstream out(outname);
  //   for(unsigned int sp=0; sp<spacepoints.size();++sp)
  //   {
  //     out << spacepoints.at(sp).x() << " "<<spacepoints.at(sp).y()<<" "<<spacepoints.at(sp).z()<<std::endl;
  //   }
  //   out.close();
  // }

  const double targetSP=10000.;
  int SPNum=spacepoints.size();
  if(SPNum<=0)
  {
    return Acts::Result<std::vector<double>>::success({-999., -999., -999., 1.0*spacepoints.size(), 0., 0.});
  }

  std::vector<float> absEtaRanges{0.3,0.6,0.9, 1.2,1.5,1.8, 2.1,2.4,2.7, 3.0};
  std::vector<float> absEtaFractions{0.05,0.05,0.06, 0.07,0.08,0.09, 0.08,0.11,0.14, 0.11};

  float maxAbsEta=3.0, minAbsEta=0.3;
  Acts::ActsScalar absEtaRange=maxAbsEta;
  Acts::ActsScalar totalFrac=0.;
  // std::cout<<" absEtaRanges.size() "<<absEtaRanges.size()<<std::endl;
  for(unsigned int r=0;r<absEtaRanges.size();++r)
  {
    Acts::ActsScalar thisRange = absEtaRanges.at(r) - (r?absEtaRanges.at(r-1):0.);
    // std::cout<<" r "<<r<<", thisRange "<<thisRange<<std::endl;
    Acts::ActsScalar addToTotalFrac = absEtaFractions.at(r);

    // std::cout<<" new total fraction "<<totalFrac<<" + "<<addToTotalFrac<<std::endl;
    // std::cout<<" which is about "<<(totalFrac+addToTotalFrac)*SPNum<<" spacepoints "<<std::endl;
    if((totalFrac+addToTotalFrac)*SPNum>targetSP)
    {
      Acts::ActsScalar needOnly = (targetSP-totalFrac*SPNum)/(SPNum*addToTotalFrac);
      // std::cout<<" needOnly "<<needOnly<<std::endl;
      absEtaRange = (r?absEtaRanges.at(r-1):0.) + needOnly*thisRange;
      // std::cout<<" absEtaRange "<<absEtaRange<<std::endl;
      break;
    }
    totalFrac += addToTotalFrac;
  }
  if(absEtaRange>maxAbsEta) absEtaRange = maxAbsEta;
  if(absEtaRange<minAbsEta) absEtaRange = minAbsEta;
  double minTheta = 2.0*std::atan(std::exp(-1.0*absEtaRange));

  double binsNumDecrease=1.;
  if(SPNum<targetSP)
  {
    binsNumDecrease = std::pow(1.35, std::log(targetSP / SPNum));
  }
  int cotThetaBins = (int)(8000./binsNumDecrease);
  
  int fillNeighbors=0, minHits=4;

  if(0.99<m_cfg.mixedEccentricity && m_cfg.mixedEccentricity<1.01)
  {
    // 4 entries SP>=1000, 3 entries SP<1000
    if(SPNum<1000) minHits=3;
  }
  if(1.99<m_cfg.mixedEccentricity && m_cfg.mixedEccentricity<2.01)
  {
    // no neighbors SP>=1000, 1 neighbor SP<1000
    if(SPNum<1000) fillNeighbors=1;
  }
  if(2.99<m_cfg.mixedEccentricity && m_cfg.mixedEccentricity<3.01)
  {
    // no neighbors SP>=1000, 1 neighbor 200<SP<1000, 3 neighbors SP<200
    if(SPNum<1000) fillNeighbors=1;
    if(SPNum<200) fillNeighbors=3;
  }
  if(3.99<m_cfg.mixedEccentricity && m_cfg.mixedEccentricity<4.01)
  {
    // no neighbors + 4 entries SP>=1000, 1 neighbor + 4 entries 200<SP<1000, 1 neighbor + 3 entries SP<200
    if(SPNum<1000) fillNeighbors=1;
    if(SPNum<200) minHits=3;
  }

  ACTS_INFO("Have m_cfg.mixedEccentricity = "<<m_cfg.mixedEccentricity<<", SPNum = "<<SPNum<<", minHits = "<<minHits<<", fillNeighbors "<<fillNeighbors);
  

  // if(SPNum*0.5 > targetSP)
  // {
  //   // for more than 40k SPs
  //   double etaRange=(targetSP/SPNum)/0.25;
  //   minTheta=2.0*std::atan(std::exp(-1.0*etaRange));

  //   if(minTheta>1.275) minTheta=1.275; // eta=0.3, for more than 333k SPs
  // }
  // else
  // {
  //   // we have little SP, make minTheta even smaller
  //   if(SPNum*0.90 < targetSP)
  //   {
  //     minTheta=0.1;  // eta=3.; for less than 22.2k SPs
  //   }
  //   else
  //   {
  //     // for 22.2k-40k SPs
  //     double etaRange=(targetSP/SPNum - 0.5)/0.4 +2.;
  //     minTheta=2.0*std::atan(std::exp(-1.0*etaRange));
  //   }
  // }


  if(m_cfg.minimalizeWRT == "hough")
  {
    double vtx_x = 0., vtx_y = 0., vtx_z;
    double step=2.0;

    // Acts::Vector2 old_vtx = {vtx_x,vtx_y};
    // ACTS_INFO("Have old vertex at "<<old_vtx[0]<<", "<<old_vtx[1]);

    auto [vtx_z1, usedSP] = findHoughPeak(spacepoints, vtx_x, vtx_y, -200, 200, minTheta, 800, cotThetaBins, fillNeighbors, minHits);
    ACTS_INFO("Have "<<spacepoints.size()<<" spacepoints, have minTheta "<<minTheta<<", actually using "<<usedSP<<" spacepoints, have "<<cotThetaBins<<" bins in cot Theta");

    auto [vtx_z2, peak_z2] = findHoughPeak(spacepoints, vtx_x, vtx_y, vtx_z1[2]-30., vtx_z1[2]+30., minTheta, 180, cotThetaBins, fillNeighbors, minHits);
    auto [vtx_z3, peak_z3] = findHoughPeak(spacepoints, vtx_x, vtx_y, vtx_z2[2]-8., vtx_z2[2]+8., minTheta, 80, cotThetaBins, fillNeighbors, minHits);
    vtx_z=vtx_z3[2];

    // ACTS_INFO("--- --- --- ---");
    // auto [vtx_z3, peak_z3] = findHoughPeak(spacepoints, 0.5,    0.,      vtx_z2[2]-8., vtx_z2[2]+8., 80);
    // ACTS_INFO("--- --- --- --- "<<vtx_z3[0]<<" "<<vtx_z3[1]<<" "<<vtx_z3[2]<<"; peak "<<peak_z3);

    for(int iter=0; iter<0; ++iter)
    {
      // ACTS_INFO("--- --- --- ---");
      auto [vtx_def, peak_def] = findHoughPeak(spacepoints, vtx_x, vtx_y, vtx_z-8., vtx_z+8., minTheta, 80, cotThetaBins, fillNeighbors, minHits);

      double new_vtx_x=vtx_x, new_vtx_y=vtx_y, new_vtx_z=vtx_z;
      double old_vtx_new_z=-999.;

      int maxGrid=(iter?3:5);
      bool haveNewVertex=false;

      for(int g=0; g<maxGrid*maxGrid; g++)
      {
        int sx=(g%maxGrid)-2;
        int sy=g/maxGrid -2;
        auto [vtx_tmp, peak_tmp] = findHoughPeak(spacepoints, 
                                                 vtx_x+sx*step, 
                                                 vtx_y+sy*step,
                                                 vtx_z-10., vtx_z+10., 
                                                 minTheta, 100, cotThetaBins, 
                                                 fillNeighbors, minHits);
        if(peak_tmp>peak_def)
        {
          new_vtx_x=vtx_x+sx*step;
          new_vtx_y=vtx_y+sy*step;
          new_vtx_z=vtx_tmp[2];
          peak_def=peak_tmp;
          haveNewVertex=true;
        }
        if(sx==0 && sy==0)
        {
          old_vtx_new_z=vtx_tmp[2];
        }
      }

      if(haveNewVertex)
      {
        vtx_x=new_vtx_x;
        vtx_y=new_vtx_y;
        vtx_z=new_vtx_z;
        // ACTS_INFO(".Have new vertex at "<<vtx_x<<", "<<vtx_y<<", "<<vtx_z<<"; new step is "<<step*0.666); 
        step*=0.666;
      }
      else
      {
        // ACTS_INFO(".Have the same vertex at "<<vtx_x<<", "<<vtx_y<<"; old z "<<vtx_z<<", new z "<<old_vtx_new_z<<"; new step is "<<step*0.333);
        vtx_z=old_vtx_new_z;
        step*=0.333;
      }

      if(step<0.1)
      {
        break;
      }
    }

    return Acts::Result<std::vector<double>>::success({vtx_x, vtx_y, vtx_z, 1.0*spacepoints.size(), usedSP, peak_z3});
  }

  return Acts::Result<std::vector<double>>::success({-999., -999., -999., 1.0*spacepoints.size(), 0., 0.});

/*
  // sort spacepoints to different phi and z slices
  Acts::SingleSeedVertexFinder<spacepoint_t>::SortedSpacepoints
      sortedSpacepoints = sortSpacepoints(spacepoints);

  // std::cout<<"spacepoints :: "<<std::endl;
  // for(std::uint32_t sp=0;sp<spacepoints.size();++sp)
  // {
  //   std::cout<<spacepoints.at(sp).x()<<" "<<
  //              spacepoints.at(sp).y()<<" "<<
  //              spacepoints.at(sp).z()<<" "<<std::endl;

  // }
  // std::cout<<"spacepoints :: "<<std::endl;

  // find triplets and fit them with plane or ray
  std::vector<Acts::SingleSeedVertexFinder<spacepoint_t>::Triplet> triplets =
      findTriplets(sortedSpacepoints);



  std::vector<std::vector<Acts::ActsScalar>> rejectVector;
  std::vector<std::vector<Acts::ActsScalar>> vtx_iter;

  std::uint32_t counter[3]={0,0,0};
  // std::cout<<"spacepoints[][3] = { "<<std::endl;
  for(std::uint32_t la=0;la<3;++la)
  {
    for(std::uint32_t phi=0;phi<m_cfg.numPhiSlices;++phi)
    {
      for(std::uint32_t z=0;z<m_cfg.numZSlices;++z)
      {
        counter[la]+=sortedSpacepoints.getSP(la,phi,z).size();
        for(std::uint32_t s=0;s<sortedSpacepoints.getSP(la,phi,z).size();++s)
        {
          // std::cout<<"{"<<
          //   sortedSpacepoints.getSP(la,phi,z).at(s).first->x()<<", "<<
          //   sortedSpacepoints.getSP(la,phi,z).at(s).first->y()<<", "<<
          //   sortedSpacepoints.getSP(la,phi,z).at(s).first->z()<<"}, ";
        }
        // std::cout<<std::endl;
      }
    }
  }
  // std::cout<<" }; "<<std::endl;
  rejectVector.push_back({1.*counter[0],1.*counter[1],1.*counter[2],1.*triplets.size()});

  // if no valid triplets found
  if (triplets.size()<m_cfg.minNumTriplets) {
    ACTS_INFO("have "<<triplets.size()<<" triplets, need at least "<<m_cfg.minNumTriplets);

    for(int i=0;i<100;++i)
    {
      vtx_iter.push_back({0.,0.,0.});
    }
    rejectVector.insert(rejectVector.begin(),vtx_iter.begin(),vtx_iter.end());

    // return {Acts::Result<Acts::Vector3>::failure(std::error_code()), rejectVector};
    return Acts::Result<Acts::Vector3>::failure(std::error_code());
  }

  // ACTS_INFO("size of 1 Triplet "<<sizeof(Triplet)<<", size of vector7 = "<<sizeof(std::vector<double>{0.,1.,2.,3.,4.,5.,6.})+sizeof(double)*7<<" size of Ray3D "<<sizeof(Acts::Ray3D));

  // ACTS_INFO("A-Size of triplets = "<<triplets.size()<<" that is "<<triplets.size()*sizeof(triplets.at(0))<<", ActsScalar "<<sizeof(Acts::ActsScalar));


  Acts::Vector3 vtx = Acts::Vector3::Zero();

  if (m_cfg.minimalizeWRT == "planes" || m_cfg.minimalizeWRT == "rays" || m_cfg.minimalizeWRT == "mixed") {
    // find a point closest to all rays fitted through the triplets
    vtx = findClosestPoint(triplets,rejectVector,vtx_iter);

    ACTS_INFO("size of SPs: "<<rejectVector.at(0).at(0)<<", "<<rejectVector.at(0).at(1)<<", "<<rejectVector.at(0).at(2)<<"; all triplets: "<<rejectVector.at(0).at(3));

    ACTS_INFO("Found mixed vertex at x = " << vtx[0]
                                     << "mm, y = " << vtx[1]
                                     << "mm, z = " << vtx[2] << "mm");

  } else {
    ACTS_ERROR("value of minimalizeWRT is "
               << m_cfg.minimalizeWRT
               << ", allowed values are \"planes\" or \"rays\" ");
  }

  // ACTS_INFO("Size of sortedSpacepoints = "<<sizeof(sortedSpacepoints.sortedSP[0][0][0])*spacepoints.size());
  // ACTS_INFO("B-Size of triplets = "<<triplets.size()<<" that is "<<triplets.size()*sizeof(triplets.at(0))*sizeof(Acts::ActsScalar)<<", ActsScalar "<<sizeof(Acts::ActsScalar));

  // FILE* file = fopen("/proc/self/status", "r");
  // char line[128];
  // while (fgets(line, 128, file) != NULL){
  //     if (strncmp(line, "VmSize:", 7) == 0){
  //         ACTS_INFO("In the end of SingleSeedVertexFinder: VmSize="<<line);
  //     }
  //     if (strncmp(line, "VmRSS:", 6) == 0){
  //       ACTS_INFO("In the end of SingleSeedVertexFinder: VmRSS="<<line);
  //     }
  // }
  // fclose(file);

  ACTS_INFO("Size of vtx_iter = "<<vtx_iter.size());
  rejectVector.insert(rejectVector.begin(),vtx_iter.begin(),vtx_iter.end());

  // return {Acts::Result<Acts::Vector3>::success(vtx), rejectVector};
  return Acts::Result<Acts::Vector3>::success(vtx);
  */
}


template <typename spacepoint_t>
std::pair<Acts::Vector3, double> Acts::SingleSeedVertexFinder<spacepoint_t>::findHoughPeak(const std::vector<spacepoint_t>& spacepoints, double vtx_x, double vtx_y, double minZ, double maxZ, double minTheta, int numZBins, int cotThetaBins, int fillNeighbors, int minHits) const {

  // int sp_counterAll=0, sp_counterSel=0, sp_counter_inBounds=0, sp_counter_outBounds=0;
  int sp_counter=0;
  // auto t1 = std::chrono::high_resolution_clock::now();

  const int sizeZ = numZBins;
  // minTheta=1.0: eta=0.6;  minTheta=0.7: eta=1.0;  minTheta=0.27: eta=2.0; 
  //   usedSP=SP*0.162    ;    usedSP=SP*0.217    ;    usedSP=SP*0.435;
  // cca. 0.2*SP per unit of |eta| will be usedSP
  // const double minTheta = 0.27;
  const double maxTheta = M_PI-minTheta;
  const double minTanTheta = std::tan(minTheta), maxTanTheta = std::tan(maxTheta);
  const double minCotTheta = 1./maxTanTheta, maxCotTheta = 1./minTanTheta;
  // const double minZ = -200., maxZ = 200.; //, maxR=400;

  const double zBinSize = (maxZ-minZ)/sizeZ;
  // const double thetaBinSize = (maxTheta-minTheta)/sizeTheta;
  const double invCotThetaBinSize = cotThetaBins/(maxCotTheta-minCotTheta);

  std::vector<int> vec_helper(cotThetaBins,0);
  std::vector<int> z_projection(sizeZ,0);
  std::vector<std::vector<int>> matrix(sizeZ,vec_helper);

  std::vector<double> vtx_z_pos;
  for(int zBin=0;zBin<sizeZ;zBin++) {
    vtx_z_pos.push_back(zBinSize * (zBin+0.5) + minZ); // vtx_x
  }

  // auto t2 = std::chrono::high_resolution_clock::now();

  for (const auto& sp : spacepoints) {
    // ++sp_counterAll;

    // if(sp.z() > maxR) continue;

    if(sp.z()>maxZ)
    {
      if(sp.r()/(sp.z()-maxZ) < minTanTheta) {
        // ACTS_INFO(sp_counterAll<<": r "<<sp.r()<<", z "<<sp.z()<<"; "<<sp.r()/(sp.z()-maxZ)<<" < "<<minTanTheta<<"; rejected by minTheta");
        continue;
      }
    }
    else if(sp.z()<minZ)
    {
      if(sp.r()/(sp.z()-minZ) > maxTanTheta) {
        // ACTS_INFO(sp_counterAll<<": r "<<sp.r()<<", z "<<sp.z()<<"; "<<sp.r()/(sp.z()-minZ)<<" > "<<maxTanTheta<<"; rejected by maxTheta");
        continue;
      }
    }
    
    // ACTS_INFO(sp_counterAll<<": r "<<sp.r()<<", z "<<sp.z()<<", accepted");

    // ++sp_counterSel;

    double sp_invr=1./std::sqrt((sp.x()-vtx_x)*(sp.x()-vtx_x) + (sp.y()-vtx_y)*(sp.y()-vtx_y));

    int zFrom = (int)(((sp.z()-1./(sp_invr*minTanTheta))-minZ)/zBinSize)+1;  
    if(zFrom<0) zFrom=0;
    int zTo   = (int)(((sp.z()-1./(sp_invr*maxTanTheta))-minZ)/zBinSize);
    if(zTo>=sizeZ) zTo=sizeZ;

    // ACTS_INFO("sp.z "<<sp.z()<<", sp_r "<<sp_r<<", minTanTheta "<<minTanTheta<<", minZ "<<minZ<<", zBinSize "<<zBinSize<<"; "<<((sp.z()-sp_r/minTanTheta)-minZ)/zBinSize<<";  zFrom "<<zFrom);
    // ACTS_INFO("sp.z "<<sp.z()<<", sp_r "<<sp_r<<", maxTanTheta "<<maxTanTheta<<", minZ "<<minZ<<", zBinSize "<<zBinSize<<"; "<<(int)(((sp.z()-sp_r/maxTanTheta)-minZ)/zBinSize)+1<<";  zTo "<<zTo);

    if(zTo>zFrom) ++sp_counter;

    for(int zBin=zFrom;zBin<zTo;zBin++) {
      // double theta = std::atan2(sp_r,(sp.z()-vtx_z));
      double cotTheta = (sp.z()-vtx_z_pos[zBin])*sp_invr;

      // int thetaBin = (int)((theta-minTheta)/thetaBinSize);
      // if(thetaBin<0 || thetaBin>=sizeTheta) {
      //   ++sp_counter_outBounds;
      //   continue;
      // }

      int cotThetaBin = (int)((cotTheta-minCotTheta)*invCotThetaBinSize);
      // if(cotThetaBin<0 || cotThetaBin>=sizeTheta) {
      //   ++sp_counter_outBounds;
      //   continue;
      // }

      // ++sp_counter_inBounds;
      for(int n=0;n<=fillNeighbors;++n)
        {
          // houghFilt[m]->Fill(vtx_z,cot+n*cotThetaBinSize[m]);
          // if(n!=0) houghFilt[m]->Fill(vtx_z,cot-n*cotThetaBinSize[m]);
          if(cotThetaBin+n < cotThetaBins) ++matrix[zBin][cotThetaBin+n];
          if(n!=0 && cotThetaBin-n >=0) ++matrix[zBin][cotThetaBin-n];
        }
    }            
  }

  // auto t3 = std::chrono::high_resolution_clock::now();

  int maxZBin=0;
  for(int zBin=0;zBin<sizeZ;zBin++) {
    // for(int thetaBin=0;thetaBin<sizeTheta;thetaBin++) {
    //   if(matrix[zBin][thetaBin] < minHits) matrix[zBin][thetaBin]=0;
    // }

    z_projection[zBin] = std::reduce(matrix[zBin].begin(), matrix[zBin].end(), 0,
    [&](int& lhs, int& rhs) { return (lhs>=minHits)*lhs + (rhs>=minHits)*rhs;} );

    // ACTS_INFO("zBin "<<zBin<<" has "<<z_projection[zBin]<<" entries; maximum at zBin "<<maxZBin<<" = "<<z_projection.at(maxZBin));
    if(z_projection[zBin] > z_projection.at(maxZBin)) maxZBin=zBin;
  }


  // auto t4 = std::chrono::high_resolution_clock::now();

  double avg=std::accumulate(z_projection.begin(), z_projection.end(), 0.)/z_projection.size();
  double sumEntries=0;
  // double onlyPeak=std::max(z_projection.at(maxZBin)-avg, 0.);
  double meanPeak=0.; //,varPeak=0., skewPeak=0.,kurtPeak=0.;

  int width=3;

  for(int zBin=std::max(maxZBin-width,0); zBin<=std::min(sizeZ-1,maxZBin+width); ++zBin) {
    sumEntries += std::max(z_projection.at(zBin)-avg, 0.);
    meanPeak   += (maxZBin-zBin)*std::max(z_projection.at(zBin)-avg, 0.);
    // varPeak    += (maxZBin-zBin)*(maxZBin-zBin)*std::max(z_projection.at(zBin)-avg, 0.);
    // skewPeak   += (maxZBin-zBin)*(maxZBin-zBin)*(maxZBin-zBin)*std::max(z_projection.at(zBin)-avg, 0.);
    // kurtPeak   += (maxZBin-zBin)*(maxZBin-zBin)*(maxZBin-zBin)*(maxZBin-zBin)*std::max(z_projection.at(zBin)-avg, 0.);
  }

 

  if(sumEntries!=0.)
  {
    meanPeak/=sumEntries;
    // varPeak/=sumEntries;
    // skewPeak/=sumEntries;
    // kurtPeak/=sumEntries;
    double real_vtx_z=(maxZ-minZ)/sizeZ * (meanPeak+maxZBin+0.5) + minZ;

    if(minZ<-199 && maxZ>199)
    {
      return {{vtx_x, vtx_y, real_vtx_z}, 1.*sp_counter};
    }

    // if(varPeak*kurtPeak != 0.) {
      return {{vtx_x, vtx_y, real_vtx_z}, sumEntries};
    // }

    // return {{vtx_x, vtx_y, real_vtx_z}, 0.};
  }

  if(minZ<-199 && maxZ>199)
  {
    return {{vtx_x, vtx_y, 0.}, 1.*sp_counter};
  }
  return {{vtx_x, vtx_y, 0.}, 0.};

  // auto t5 = std::chrono::high_resolution_clock::now();

  // ACTS_INFO("Initialization done in " << (t2 - t1).count() / 1e6 << " ms");
  // ACTS_INFO("Matrix filled  in "      << (t3 - t2).count() / 1e6 << " ms; actual SPs "<<sp_counterAll<<" / "<<sp_counterSel<<"; in bouns "<<sp_counter_inBounds<<" / out "<<sp_counter_outBounds);
  // ACTS_INFO("Nullification done in "  << (t4 - t3).count() / 1e6 << " ms");
  // ACTS_INFO("Final peak done in "     << (t5 - t4).count() / 1e6 << " ms"<<"; at "<<vtx_x<<" "<<vtx_y<<" "<<real_vtx_z<<"; peak "<<sumEntries<<", varPeak "<<varPeak<<", skewPeak "<<skewPeak<<", kurtPeak "<<kurtPeak<<"; sumEntries/(varPeak*kurtPeak) "<<sumEntries/(varPeak*kurtPeak)<<", peak/varPeak "<<sumEntries/varPeak);

  // return {{vtx_x, vtx_y, real_vtx_z}, sumEntries/(varPeak*kurtPeak)};
}


template <typename spacepoint_t>
typename Acts::SingleSeedVertexFinder<spacepoint_t>::SortedSpacepoints
Acts::SingleSeedVertexFinder<spacepoint_t>::sortSpacepoints(
    const std::vector<spacepoint_t>& spacepoints) const {
  Acts::SingleSeedVertexFinder<spacepoint_t>::SortedSpacepoints
      sortedSpacepoints(m_cfg.numPhiSlices, m_cfg.numZSlices);

  for (const auto& sp : spacepoints) {
    // phi will be saved for later
    Acts::ActsScalar phi = detail::radian_pos(std::atan2(sp.y(), sp.x()));
    std::uint32_t phislice =
        (std::uint32_t)(phi / (2 * M_PI) * m_cfg.numPhiSlices);
    if (phislice >= m_cfg.numPhiSlices) {
      phislice = 0;
    }

    if (std::abs(sp.z()) >= m_cfg.maxAbsZ) {
      continue;
    }
    std::uint32_t zslice = (std::uint32_t)(
        (sp.z() + m_cfg.maxAbsZ) / (2 * m_cfg.maxAbsZ) * m_cfg.numZSlices);

    // input spacepoint is sorted into one of the subsets
    if (sp.r() < m_cfg.rMinMiddle) {
      if (m_cfg.rMinNear < sp.r() && sp.r() < m_cfg.rMaxNear) {
        if (std::fmod(m_cfg.useFracPhiSlices * phislice, 1.0) >=
            m_cfg.useFracPhiSlices) {
          continue;
        }
        sortedSpacepoints.addSP(0, phislice, zslice)
            .emplace_back((spacepoint_t const*)&sp, phi);
      }
    } else if (sp.r() < m_cfg.rMinFar) {
      if (sp.r() < m_cfg.rMaxMiddle) {
        if (std::fmod(m_cfg.useFracZSlices * zslice, 1.0) >=
            m_cfg.useFracZSlices) {
          continue;
        }
        sortedSpacepoints.addSP(1, phislice, zslice)
            .emplace_back((spacepoint_t const*)&sp, phi);
      }
    } else if (sp.r() < m_cfg.rMaxFar) {
      sortedSpacepoints.addSP(2, phislice, zslice)
          .emplace_back((spacepoint_t const*)&sp, phi);
    }
  }

  std::uint32_t counter[3]={0,0,0};
  for(std::uint32_t la=0;la<3;++la)
  {
    for(std::uint32_t phi=0;phi<m_cfg.numPhiSlices;++phi)
    {
      for(std::uint32_t z=0;z<m_cfg.numZSlices;++z)
      {
        counter[la]+=sortedSpacepoints.getSP(la,phi,z).size();
      }
    }
  }
  ACTS_INFO("have sorted spacepoints: "<<counter[0]<<", "<<counter[1]<<", "<<counter[2]);

  return sortedSpacepoints;
}

template <typename spacepoint_t>
std::vector<typename Acts::SingleSeedVertexFinder<spacepoint_t>::Triplet>
Acts::SingleSeedVertexFinder<spacepoint_t>::findTriplets(
    const Acts::SingleSeedVertexFinder<spacepoint_t>::SortedSpacepoints&
        sortedSpacepoints) const {
  std::vector<Acts::SingleSeedVertexFinder<spacepoint_t>::Triplet> triplets;

  std::uint32_t phiStep =
      (std::uint32_t)(m_cfg.maxPhideviation / (2 * M_PI / m_cfg.numPhiSlices)) +
      1;

  std::uint32_t considered_triplets = 0;
  std::uint32_t rej0=0, rej1 = 0, rej2 = 0, rej3 = 0, rej4 = 0;
  std::uint32_t good_triplets = 0;

  // calculate limits for middle spacepoints
  Acts::Vector2 vecA{-m_cfg.maxAbsZ + m_cfg.maxZPosition, m_cfg.rMinFar};
  vecA /= 2.;
  Acts::Vector2 vecB = {vecA[1], -vecA[0]};
  vecB /= std::tan(m_cfg.maxXYZdeviation);
  Acts::Vector2 posR = Acts::Vector2(-m_cfg.maxZPosition, 0.) + vecA + vecB;
  Acts::ActsScalar R = vecA.norm() / std::sin(m_cfg.maxXYZdeviation);
  Acts::ActsScalar constB = -2. * posR[0];
  Acts::ActsScalar constC =
      posR[0] * posR[0] +
      (posR[1] - m_cfg.rMaxNear) * (posR[1] - m_cfg.rMaxNear) - R * R;
  Acts::ActsScalar maxZMiddle =
      -1. * (-constB - sqrt(constB * constB - 4. * constC)) / 2.;
  if (maxZMiddle <= 0) {
    ACTS_WARNING(
        "maximum position of middle spacepoints is not positive, maxZMiddle = "
        << maxZMiddle << ", check your config; setting maxZMiddle to "
        << m_cfg.maxAbsZ);
    maxZMiddle = m_cfg.maxAbsZ;
  }

  // save some constant values for later
  Acts::ActsScalar rNearRatio[2] = {m_cfg.rMinNear / m_cfg.rMaxMiddle,
                                    m_cfg.rMaxNear / m_cfg.rMinMiddle};
  Acts::ActsScalar rMiddle[2] = {m_cfg.rMaxMiddle, m_cfg.rMinMiddle};
  Acts::ActsScalar rFarDelta[2] = {
      m_cfg.rMaxFar - m_cfg.rMinMiddle,
      m_cfg.rMinFar - m_cfg.rMaxMiddle,
  };
  Acts::ActsScalar zBinLength = 2. * m_cfg.maxAbsZ / m_cfg.numZSlices;

  // limits in terms of slice numbers
  std::uint32_t limitMiddleSliceFrom =
      (std::uint32_t)((-maxZMiddle + m_cfg.maxAbsZ) / zBinLength);
  std::uint32_t limitMiddleSliceTo =
      (std::uint32_t)((maxZMiddle + m_cfg.maxAbsZ) / zBinLength + 1);
  std::uint32_t limitAbsZSliceFrom = (std::uint32_t)(
      (-m_cfg.maxZPosition + m_cfg.maxAbsZ) / zBinLength + 0.01);
  std::uint32_t limitAbsZSliceTo =
      (std::uint32_t)((m_cfg.maxZPosition + m_cfg.maxAbsZ) / zBinLength + 1.01);

  for (std::uint32_t middleZ = limitMiddleSliceFrom;
       middleZ < limitMiddleSliceTo; ++middleZ) {
    // skip slices that are empty anyway
    if (std::fmod(m_cfg.useFracZSlices * middleZ, 1.0) >=
        m_cfg.useFracZSlices) {
      continue;
    }

    // calculate limits for near spacepoints, assuming the middle spacepoints
    // are within some boundaries
    bool isLessFrom = (middleZ <= limitAbsZSliceFrom);
    Acts::ActsScalar deltaZfrom =
        (middleZ - limitAbsZSliceFrom - 1) * zBinLength;
    Acts::ActsScalar angleZfrom =
        std::atan2(rMiddle[isLessFrom], deltaZfrom) + m_cfg.maxXYZdeviation;
    std::uint32_t nearZFrom = 0;
    if (angleZfrom < M_PI) {
      Acts::ActsScalar new_deltaZfrom =
          rMiddle[isLessFrom] / std::tan(angleZfrom) / zBinLength;
      nearZFrom = (std::uint32_t)std::max(
          new_deltaZfrom * rNearRatio[isLessFrom] + limitAbsZSliceFrom, 0.);
    }

    bool isLessTo = (middleZ < limitAbsZSliceTo);
    Acts::ActsScalar deltaZto = (middleZ - limitAbsZSliceTo + 1) * zBinLength;
    Acts::ActsScalar angleZto =
        std::atan2(rMiddle[!isLessTo], deltaZto) - m_cfg.maxXYZdeviation;
    std::uint32_t nearZTo = m_cfg.numZSlices;
    if (angleZto > 0) {
      Acts::ActsScalar new_deltaZto =
          rMiddle[!isLessTo] / std::tan(angleZto) / zBinLength;
      nearZTo = (std::uint32_t)std::max(
          new_deltaZto * rNearRatio[!isLessTo] + limitAbsZSliceTo, 0.);
      if (nearZTo > m_cfg.numZSlices) {
        nearZTo = m_cfg.numZSlices;
      }
    }

    for (std::uint32_t nearZ = nearZFrom; nearZ < nearZTo; ++nearZ) {
      // calculate limits for far spacepoits, assuming middle and near
      // spacepoits are within some boundaries
      bool isMiddleLess = (middleZ <= nearZ);

      Acts::ActsScalar delta2Zfrom = (middleZ - nearZ - 1) * zBinLength;
      Acts::ActsScalar angle2Zfrom =
          std::atan2(rFarDelta[isMiddleLess], delta2Zfrom) +
          m_cfg.maxXYZdeviation;
      std::uint32_t farZFrom = 0;
      if (angle2Zfrom < M_PI) {
        farZFrom = (std::uint32_t)std::max(
            (rFarDelta[isMiddleLess] / std::tan(angle2Zfrom) / zBinLength) +
                middleZ,
            0.);
        if (farZFrom >= m_cfg.numZSlices) {
          continue;
        }
      }

      isMiddleLess = (middleZ < nearZ);
      Acts::ActsScalar delta2Zto = (middleZ - nearZ + 1) * zBinLength;
      Acts::ActsScalar angle2Zto =
          std::atan2(rFarDelta[!isMiddleLess], delta2Zto) -
          m_cfg.maxXYZdeviation;
      std::uint32_t farZTo = m_cfg.numZSlices;
      if (angle2Zto > 0) {
        farZTo = (std::uint32_t)std::max(
            (rFarDelta[!isMiddleLess] / std::tan(angle2Zto) / zBinLength) +
                middleZ + 1,
            0.);
        if (farZTo > m_cfg.numZSlices) {
          farZTo = m_cfg.numZSlices;
        } else if (farZTo == 0) {
          continue;
        }
      }

      for (std::uint32_t farZ = farZFrom; farZ < farZTo; farZ++) {
        // loop over near phi slices
        for (std::uint32_t nearPhi = 0; nearPhi < m_cfg.numPhiSlices;
             ++nearPhi) {
          // skip slices that are empty anyway
          if (std::fmod(m_cfg.useFracPhiSlices * nearPhi, 1.0) >=
              m_cfg.useFracPhiSlices) {
            continue;
          }

          // loop over some middle phi slices
          for (std::uint32_t middlePhi_h =
                   m_cfg.numPhiSlices + nearPhi - phiStep;
               middlePhi_h <= m_cfg.numPhiSlices + nearPhi + phiStep;
               ++middlePhi_h) {
            std::uint32_t middlePhi = middlePhi_h % m_cfg.numPhiSlices;

            // loop over some far phi slices
            for (std::uint32_t farPhi_h =
                     m_cfg.numPhiSlices + middlePhi - phiStep;
                 farPhi_h <= m_cfg.numPhiSlices + middlePhi + phiStep;
                 ++farPhi_h) {
              std::uint32_t farPhi = farPhi_h % m_cfg.numPhiSlices;

              // for all near spacepoints in this slice
              for (const auto& nearSP :
                   sortedSpacepoints.getSP(0, nearPhi, nearZ)) {
                Acts::ActsScalar phiA = nearSP.second;

                // for all middle spacepoints in this slice
                for (const auto& middleSP :
                     sortedSpacepoints.getSP(1, middlePhi, middleZ)) {
                  Acts::ActsScalar phiB = middleSP.second;
                  Acts::ActsScalar deltaPhiAB =
                      detail::difference_periodic(phiA, phiB, 2 * M_PI);
                  if (std::abs(deltaPhiAB) > m_cfg.maxPhideviation) {
                    continue;
                  }

                  // for all far spacepoints in this slice
                  for (const auto& farSP :
                       sortedSpacepoints.getSP(2, farPhi, farZ)) {
                    Acts::ActsScalar phiC = farSP.second;
                    Acts::ActsScalar deltaPhiBC =
                        detail::difference_periodic(phiB, phiC, 2 * M_PI);
                    if (std::abs(deltaPhiBC) > m_cfg.maxPhideviation) {
                      continue;
                    }

                    auto fittedTriplet = tripletValidationAndFit(
                        *nearSP.first, *middleSP.first, *farSP.first);

                    ++considered_triplets;

                    if (fittedTriplet.getDistance()>-1.5) {
                      triplets.push_back(fittedTriplet);
                      ++good_triplets;
                    }
                    else
                    {
                      ++rej0;
                      if(fittedTriplet.getDistance()>-2.15) ++rej1;
                      if(fittedTriplet.getDistance()>-2.25) ++rej2;
                      if(fittedTriplet.getDistance()>-2.35) ++rej3;
                      if(fittedTriplet.getDistance()>-2.45) ++rej4;
                    }
                  }  // loop over far spacepoints
                }    // loop over middle spacepoints
              }      // loop over near spacepoints
            }        // loop over far phi slices
          }          // loop over middle phi slices
        }            // loop over near phi slices
      }              // loop over far Z slices
    }                // loop over near Z slices
  }                  // loop over middle Z slices

  ACTS_INFO("considered "<<considered_triplets<<", good "<<good_triplets<<", i.e. bad "<<considered_triplets-good_triplets);
  ACTS_INFO("rejected  total - "<<rej0<<"; rejected by step 1 - "<<rej1<<", by step 2 - "<<rej2<<", by step 3 - "<<rej3<<", by step 4 - "<<rej4);
  return triplets;
}

template <typename spacepoint_t>
typename Acts::SingleSeedVertexFinder<spacepoint_t>::Triplet Acts::SingleSeedVertexFinder<spacepoint_t>::tripletValidationAndFit(
    const spacepoint_t& a, const spacepoint_t& b, const spacepoint_t& c) const {
  // slope for near+middle spacepoints
  Acts::ActsScalar alpha1 =
      std::atan2(a.y() - b.y(), a.x() - b.x());
  // slope for middle+far spacepoints
  Acts::ActsScalar alpha2 =
      std::atan2(b.y() - c.y(), b.x() - c.x());
  // these two slopes shouldn't be too different
  Acts::ActsScalar deltaAlpha =
      detail::difference_periodic(alpha1, alpha2, 2 * M_PI);
  if (std::abs(deltaAlpha) > m_cfg.maxXYdeviation) {
    return Triplet(-2.1, a,b,c);
  }

  // near-middle ray
  Acts::Vector3 ab{a.x() - b.x(), a.y() - b.y(),
                   a.z() - b.z()};
  // middle-far ray
  Acts::Vector3 bc{b.x() - c.x(), b.y() - c.y(),
                   b.z() - c.z()};
  // dot product of these two
  Acts::ActsScalar cosTheta = (ab.dot(bc)) / (ab.norm() * bc.norm());
  Acts::ActsScalar theta = std::acos(cosTheta);
  if (theta > m_cfg.maxXYZdeviation) {
    return Triplet(-2.2, a,b,c);
  }

  // reject the ray if it doesn't come close to the z-axis
  auto triplet = fitTriplet(a, b, c);
  const Acts::Vector3& startPoint = triplet.getStartPoint();
  const Acts::Vector3& direction = triplet.getDirection();
  // ACTS_INFO("fitted triplet: "<<startPoint[0]<<" "<<startPoint[1]<<" "<<startPoint[2]<<" .. "<<direction[0]<<" "<<direction[1]<<" "<<direction[2]);
  // norm to z-axis and to the ray
  Acts::Vector3 norm{-1. * direction[1], 1. * direction[0], 0};
  Acts::ActsScalar norm_size = norm.norm();

  Acts::ActsScalar tanTheta = norm_size / direction[2];
  if (std::abs(tanTheta) < std::tan(m_cfg.minTheta)) {
    return Triplet(-2.3,a,b,c);
  }

  // nearest distance from the ray to z-axis
  Acts::ActsScalar dist = std::abs(startPoint.dot(norm)) / norm_size;
  if (dist > m_cfg.maxRPosition) {
    return Triplet(-2.4, a,b,c);
  }

  // z coordinate of the nearest distance from the ray to z-axis
  Acts::ActsScalar zDist =
      direction.cross(norm).dot(startPoint) / (norm_size * norm_size);
  if (std::abs(zDist) > m_cfg.maxZPosition) {
    return Triplet(-2.5, a,b,c);
  }

  return triplet;
}

template <typename spacepoint_t>
typename Acts::SingleSeedVertexFinder<spacepoint_t>::Triplet
Acts::SingleSeedVertexFinder<spacepoint_t>::fitTriplet(
    const spacepoint_t& spA, const spacepoint_t& spB, const spacepoint_t& spC) const {
  // always fit ray
  auto ray = makeRayFromTriplet(spA,spB,spC);

  std::pair<Acts::Vector3, Acts::ActsScalar> plane{{0.,0.,1.},0.};
  if(m_effectEccSq!=0.)
  {
    // fit plane only if minimalizeWRT is not "rays"
    plane = makePlaneFromTriplet(spA,spB,spC);
  }

  return Triplet(plane.first, plane.second, ray.first, ray.second, spA,spB,spC);
}

template <typename spacepoint_t>
std::pair<Acts::Vector3, Acts::ActsScalar>
Acts::SingleSeedVertexFinder<spacepoint_t>::makePlaneFromTriplet(
    const spacepoint_t& spA, const spacepoint_t& spB, const spacepoint_t& spC) const {
  Acts::Vector3 a{spA.x(), spA.y(), spA.z()};
  Acts::Vector3 b{spB.x(), spB.y(), spB.z()};
  Acts::Vector3 c{spC.x(), spC.y(), spC.z()};

  Acts::Vector3 ba = b - a, ca = c - a;

  // vector (alpha,beta,gamma) normalized to unity for convenience
  Acts::Vector3 abg = ba.cross(ca).normalized();
  Acts::ActsScalar delta = -1. * abg.dot(a);

  // plane (alpha*x + beta*y + gamma*z + delta = 0)
  return {abg, delta};
}

template <typename spacepoint_t>
std::pair<Acts::Vector3, Acts::Vector3> Acts::SingleSeedVertexFinder<spacepoint_t>::makeRayFromTriplet(
    const spacepoint_t& a, const spacepoint_t& b, const spacepoint_t& c) const {
  Acts::SquareMatrix3 mat;
  mat.row(0) = Acts::Vector3(a.x(), a.y(), a.z());
  mat.row(1) = Acts::Vector3(b.x(), b.y(), b.z());
  mat.row(2) = Acts::Vector3(c.x(), c.y(), c.z());

  Acts::Vector3 mean = mat.colwise().mean();
  Acts::SquareMatrix3 cov = (mat.rowwise() - mean.transpose()).transpose() *
                         (mat.rowwise() - mean.transpose()) / 3.;

  // "cov" is self-adjoint matrix
  Eigen::SelfAdjointEigenSolver<Acts::SquareMatrix3> saes(cov);
  // eigenvalues are sorted in increasing order
  Acts::Vector3 eivec = saes.eigenvectors().col(2);

  return {mean, eivec};
}


template <typename spacepoint_t>
Acts::Vector3
Acts::SingleSeedVertexFinder<spacepoint_t>::findClosestPoint(
    std::vector<typename Acts::SingleSeedVertexFinder<spacepoint_t>::Triplet>& allTriples,
    std::vector<std::vector<Acts::ActsScalar>>& rejectVector,
    std::vector<std::vector<Acts::ActsScalar>>& vtx_iter) const {
  // 1. define function f = sum over all triplets [distance from an unknown
  // point
  //    (x_0,y_0,z_0) to the plane defined by the triplet]
  // 2. find minimum of "f" by partial derivations over x_0, y_0, and z_0
  // 3. each derivation has parts linearly depending on x_0, y_0, and z_0
  //    (will fill A[deriv][3]) or to nothing (will fill B[deriv])
  // 4. solve A*(x_0,y_0,z_0) = B

  Acts::Vector3 vtx = Acts::Vector3::Zero();
  Acts::Vector3 vtxPrev{m_cfg.rMaxFar, m_cfg.rMaxFar, m_cfg.maxAbsZ};
  // rejectVector.clear();
  vtx_iter.clear();

  // ACTS_INFO("A-size of Acts::Vector3 "<<sizeof(Acts::Vector3)<<", "<<sizeof(std::pair<Acts::Vector3, Acts::ActsScalar>)<<"; tripletsWithRays pairs "<<sizeof(std::pair<std::pair<Acts::Vector3, Acts::ActsScalar>, Acts::ActsScalar>)<<", size "<<allTriples.size());

  ACTS_INFO("A-m_effectEccSq "<<m_effectEccSq);

  // elements of the linear equations to solve
  Acts::SquareMatrix3 A = Acts::SquareMatrix3::Zero() + (1.-m_effectEccSq)*Acts::SquareMatrix3::Identity() * 2. * allTriples.size();
  Acts::Vector3 B = Acts::Vector3::Zero();
  for (const auto& triplet : allTriples) {
    const Acts::Vector3& abg = triplet.getPlaneABG();
    const Acts::ActsScalar& delta = triplet.getPlaneDelta();
    const Acts::Vector3& startPoint=triplet.getStartPoint();
    const Acts::Vector3& direction=triplet.getDirection();

    A += m_effectEccSq * 2. * (abg * abg.transpose()) - (1-m_effectEccSq) * 2. * (direction * direction.transpose());
    B -= m_effectEccSq * 2. * delta * abg + (1-m_effectEccSq) * (2. * direction * (direction.dot(startPoint)) - 2. * startPoint);
  }

  for (std::uint32_t iter = 0; iter <= m_cfg.maxIterations; iter++) {
    // new vertex position
    vtx = A.lu().solve(B);
    vtx_iter.push_back({vtx[0],vtx[1],vtx[2]});

    Acts::Vector3 vtxDiff = vtx - vtxPrev;

    ACTS_INFO(iter<<": vtx = "<<vtx[0]<<", "<<vtx[1]<<", "<<vtx[2]<<"; vtxDiff = "<<vtxDiff[0]<<", "<<vtxDiff[1]<<", "<<vtxDiff[2]<<"; there are "<<allTriples.size()<<" triplets");

    if (iter>=m_cfg.minIterations && vtxDiff.norm() < m_cfg.minVtxShift) {
      // difference between the new vertex and the old vertex is not so large
      break;
    }

    if (iter != m_cfg.maxIterations) {
      // is not the last iteration
      vtxPrev = vtx;
      // int cnt=0;
      for (auto& triplet : allTriples) {
        const Acts::Vector3& abg = triplet.getPlaneABG();
        const Acts::ActsScalar& delta = triplet.getPlaneDelta();
        const Acts::Vector3& startPoint=triplet.getStartPoint();
        const Acts::Vector3& direction=triplet.getDirection();

        // distance from plane ^2  +  distance from ray ^2
        // ACTS_INFO(iter<<" "<<cnt<<" : abg = "<<abg[0]<<", "<<abg[1]<<", "<<abg[2]<<"; delta = "<<delta<<"; startPoint "<<startPoint[0]<<", "<<startPoint[1]<<", "<<startPoint[2]<<", direction = "<<direction[0]<<", "<<direction[1]<<", "<<direction[2]);
        // ACTS_INFO(iter<<" "<<cnt<<" : m_effectEccSq = "<<m_effectEccSq<<", dist sq to plane "<<std::abs(abg.dot(vtx) + delta) * std::abs(abg.dot(vtx) + delta) <<", dist sq to ray "<<(vtx - startPoint).cross(direction).squaredNorm());


        const Acts::ActsScalar distanceSq = m_effectEccSq * std::abs(abg.dot(vtx) + delta) * std::abs(abg.dot(vtx) + delta) + (1.-m_effectEccSq) * (vtx - startPoint).cross(direction).squaredNorm();
        triplet.setDistance(distanceSq);
        // ACTS_INFO(iter<<" "<<cnt<<" : distance = "<<std::sqrt(distanceSq)<<", distance for rays = "<<std::sqrt((vtx - startPoint).cross(direction).squaredNorm()));
        // ++cnt;
      }

      std::sort(allTriples.begin(), allTriples.end(),
                [](const auto& lhs, const auto& rhs) {
                  return lhs.getDistance() < rhs.getDistance();
                });

      std::uint32_t threshold = (std::uint32_t)(allTriples.size() *
                                                (1. - m_cfg.removeFraction));

      for (std::uint32_t tr = threshold + 1; tr < allTriples.size();
           ++tr) {
        const Acts::Vector3& abg = allTriples[tr].getPlaneABG();
        const Acts::ActsScalar& delta = allTriples[tr].getPlaneDelta();
        const Acts::Vector3& startPoint = allTriples[tr].getStartPoint();
        const Acts::Vector3& direction  = allTriples[tr].getDirection();

        rejectVector.push_back({iter*1.+m_effectEccSq*0.1,abg[0],abg[1],abg[2],delta,
                                startPoint[0],startPoint[1],startPoint[2],
                                direction[0],direction[1],direction[2]});

        // remove this triplet from A and B
        A -= m_effectEccSq * 2. * (abg * abg.transpose()) + (1-m_effectEccSq) * (Acts::SquareMatrix3::Identity() * 2. -  2. * (direction * direction.transpose()));
        B += m_effectEccSq * 2. * delta * abg + (1-m_effectEccSq) * (2. * direction * (direction.dot(startPoint)) - 2. * startPoint);
      }

      // remove all excessive triplets
      if(threshold + 1 < allTriples.size()) {
        ACTS_INFO("removed "<<allTriples.size()-threshold-1<<" triplets with removeFraction");
        allTriples.erase(allTriples.begin()+threshold+1,allTriples.end());
      }

      // std::uint32_t remove_counter=0;
      // std::uint32_t good_threshold = (std::uint32_t)(tripletsWithRays.size()*0.40);
      // std::uint32_t bad_threshold = (std::uint32_t)(tripletsWithRays.size()*0.80);
      // for (std::uint32_t tr = 0; tr < good_threshold; ++tr) {

      //   auto& good_triplet = tripletsWithRays[tr];

      //   for (std::uint32_t tr2 = std::max(tr+1,bad_threshold); tr2 < tripletsWithRays.size(); ++tr2) {
      //     auto& bad_triplet = tripletsWithRays[tr2];

      //     if(good_triplet.nearSP==bad_triplet.nearSP || 
      //        good_triplet.middleSP==bad_triplet.middleSP || 
      //        good_triplet.farSP==bad_triplet.farSP) {
      //       //ACTS_INFO("good triplet tr "<<tr<<" has overlap with bad triplet tr2 "<<tr2);

      //       const Acts::Vector3& abg = tripletsWithPlanes[tr2].getPlaneABG();
      //       const Acts::ActsScalar& delta = tripletsWithPlanes[tr2].getPlaneDelta();

      //       // remove this triplet from A and B
      //       A -= 2. * (abg * abg.transpose());
      //       B += 2. * delta * abg;

      //       tripletsWithPlanes.erase(tripletsWithPlanes.begin()+tr2,tripletsWithPlanes.begin()+tr2+1);

      //       --tr2;
      //       ++remove_counter;
      //     }
      //   }
      // }
      // ACTS_INFO("removed "<<remove_counter<<" triplets with good-bad triplets / 1");

      // remove_counter=0;
      // good_threshold = (std::uint32_t)(tripletsWithPlanes.size()*0.90);
      // bad_threshold = (std::uint32_t)(tripletsWithPlanes.size()*0.60);
      // for (std::uint32_t tr = 0; tr < good_threshold; ++tr) {

      //   auto& good_triplet = tripletsWithPlanes[tr];

      //   for (std::uint32_t tr2 = std::max(tr+1,bad_threshold); tr2 < tripletsWithPlanes.size(); ++tr2) {
      //     auto& bad_triplet = tripletsWithPlanes[tr2];

      //     if((good_triplet.nearSP==bad_triplet.nearSP && 
      //         (good_triplet.middleSP==bad_triplet.middleSP || good_triplet.farSP==bad_triplet.farSP)) ||
      //         (good_triplet.middleSP==bad_triplet.middleSP && good_triplet.farSP==bad_triplet.farSP)) {
      //       //ACTS_INFO("good triplet tr "<<tr<<" has overlap with bad triplet tr2 "<<tr2);

      //       const Acts::Vector3& abg = tripletsWithPlanes[tr2].getPlaneABG();
      //       const Acts::ActsScalar& delta = tripletsWithPlanes[tr2].getPlaneDelta();

      //       // remove this triplet from A and B
      //       A -= 2. * (abg * abg.transpose());
      //       B += 2. * delta * abg;

      //       tripletsWithPlanes.erase(tripletsWithPlanes.begin()+tr2,tripletsWithPlanes.begin()+tr2+1);

      //       --tr2;
      //       ++remove_counter;
      //     }
      //   }
      // }
      // ACTS_INFO("removed "<<remove_counter<<" triplets with good-bad triplets / 2");


    }
  }

  for (std::uint32_t iter = vtx_iter.size(); iter < 99; iter++) {
    vtx_iter.push_back({0.,0.,0.});
  }
  vtx_iter.push_back({vtx[0],vtx[1],vtx[2]});

  ACTS_INFO("after fills, vtx_iter.size = "<<vtx_iter.size());

  for (std::uint32_t tr = 0; tr < allTriples.size();++tr)
  {
    const Acts::Vector3& abg = allTriples[tr].getPlaneABG();
    const Acts::ActsScalar& delta = allTriples[tr].getPlaneDelta();
    const Acts::Vector3& startPoint = allTriples[tr].getStartPoint();
    const Acts::Vector3& direction  = allTriples[tr].getDirection();

    rejectVector.push_back({99.+m_effectEccSq*0.1,abg[0],abg[1],abg[2],delta,
                            startPoint[0],startPoint[1],startPoint[2],
                            direction[0],direction[1],direction[2]});
  }

  ACTS_INFO("B-size of Acts::Vector3 "<<sizeof(Acts::Vector3)<<", "<<sizeof(std::pair<Acts::Vector3, Acts::ActsScalar>)<<"; allTriples pairs "<<sizeof(std::pair<std::pair<Acts::Vector3, Acts::ActsScalar>, Acts::ActsScalar>)<<", size "<<allTriples.size());

  return vtx;
}

////////////////////////////////

/*
template <typename spacepoint_t>
Acts::Vector3
Acts::SingleSeedVertexFinder<spacepoint_t>::findClosestPointFromPlanes(
    std::vector<typename Acts::SingleSeedVertexFinder<spacepoint_t>::Triplet>& tripletsWithPlanes) const {
  // 1. define function f = sum over all triplets [distance from an unknown
  // point
  //    (x_0,y_0,z_0) to the plane defined by the triplet]
  // 2. find minimum of "f" by partial derivations over x_0, y_0, and z_0
  // 3. each derivation has parts linearly depending on x_0, y_0, and z_0
  //    (will fill A[deriv][3]) or to nothing (will fill B[deriv])
  // 4. solve A*(x_0,y_0,z_0) = B
  Acts::Vector3 vtx = Acts::Vector3::Zero();
  Acts::Vector3 vtxPrev{m_cfg.rMaxFar, m_cfg.rMaxFar, m_cfg.maxAbsZ};
  ACTS_INFO("A-size of Acts::Vector3 "<<sizeof(Acts::Vector3)<<", "<<sizeof(std::pair<Acts::Vector3, Acts::ActsScalar>)<<"; tripletsWithRays pairs "<<sizeof(std::pair<std::pair<Acts::Vector3, Acts::ActsScalar>, Acts::ActsScalar>)<<", size "<<tripletsWithPlanes.size());
  // elements of the linear equations to solve
  Acts::SquareMatrix3 A = Acts::SquareMatrix3::Zero();
  Acts::Vector3 B = Acts::Vector3::Zero();
  for (const auto& triplet : tripletsWithPlanes) {
    const Acts::Vector3& abg = triplet.getPlaneABG();
    const Acts::ActsScalar& delta = triplet.getPlaneDelta();
    A += 2. * (abg * abg.transpose());
    B -= 2. * delta * abg;
  }
  for (std::uint32_t iter = 0; iter <= m_cfg.maxIterations; iter++) {
    // new vertex position
    vtx = A.lu().solve(B);
    Acts::Vector3 vtxDiff = vtx - vtxPrev;
    ACTS_INFO(iter<<": vtx = "<<vtx[0]<<", "<<vtx[1]<<", "<<vtx[2]<<"; vtxDiff = "<<vtxDiff[0]<<", "<<vtxDiff[1]<<", "<<vtxDiff[2]<<"; there are "<<tripletsWithPlanes.size()<<" triplets");
    if (vtxDiff.norm() < m_cfg.minVtxShift) {
      // difference between the new vertex and the old vertex is not so large
      break;
    }
    if (iter != m_cfg.maxIterations) {
      // is not the last iteration
      vtxPrev = vtx;
      for (auto& triplet : tripletsWithPlanes) {
        const Acts::Vector3& abg = triplet.getPlaneABG();
        const Acts::ActsScalar& delta = triplet.getPlaneDelta();
        const Acts::ActsScalar distance = std::abs(abg.dot(vtx) + delta);
        triplet.setDistance(distance);
      }
      std::sort(tripletsWithPlanes.begin(), tripletsWithPlanes.end(),
                [](const auto& lhs, const auto& rhs) {
                  return lhs.getDistance() < rhs.getDistance();
                });
      std::uint32_t threshold = (std::uint32_t)(tripletsWithPlanes.size() *
                                                (1. - m_cfg.removeFraction));
      for (std::uint32_t tr = threshold + 1; tr < tripletsWithPlanes.size();
           ++tr) {
        const Acts::Vector3& abg = tripletsWithPlanes[tr].getPlaneABG();
        const Acts::ActsScalar& delta = tripletsWithPlanes[tr].getPlaneDelta();
        // remove this triplet from A and B
        A -= 2. * (abg * abg.transpose());
        B += 2. * delta * abg;
      }
      // remove all excessive triplets
      if(threshold + 1 < tripletsWithPlanes.size()) {
        ACTS_INFO("removed "<<tripletsWithPlanes.size()-threshold-1<<" triplets with removeFraction");
        tripletsWithPlanes.erase(tripletsWithPlanes.begin()+threshold+1,tripletsWithPlanes.end());
      }
*/
/*
      std::uint32_t remove_counter=0;
      std::uint32_t good_threshold = (std::uint32_t)(tripletsWithPlanes.size()*0.40);
      std::uint32_t bad_threshold = (std::uint32_t)(tripletsWithPlanes.size()*0.80);
      for (std::uint32_t tr = 0; tr < good_threshold; ++tr) {
        auto& good_triplet = tripletsWithPlanes[tr];
        for (std::uint32_t tr2 = std::max(tr+1,bad_threshold); tr2 < tripletsWithPlanes.size(); ++tr2) {
          auto& bad_triplet = tripletsWithPlanes[tr2];
          if(good_triplet.nearSP==bad_triplet.nearSP || 
             good_triplet.middleSP==bad_triplet.middleSP || 
             good_triplet.farSP==bad_triplet.farSP) {
            //ACTS_INFO("good triplet tr "<<tr<<" has overlap with bad triplet tr2 "<<tr2);
            const Acts::Vector3& abg = tripletsWithPlanes[tr2].getPlaneABG();
            const Acts::ActsScalar& delta = tripletsWithPlanes[tr2].getPlaneDelta();
            // remove this triplet from A and B
            A -= 2. * (abg * abg.transpose());
            B += 2. * delta * abg;
            tripletsWithPlanes.erase(tripletsWithPlanes.begin()+tr2,tripletsWithPlanes.begin()+tr2+1);
            --tr2;
            ++remove_counter;
          }
        }
      }
      ACTS_INFO("removed "<<remove_counter<<" triplets with good-bad triplets / 1");
      remove_counter=0;
      good_threshold = (std::uint32_t)(tripletsWithPlanes.size()*0.90);
      bad_threshold = (std::uint32_t)(tripletsWithPlanes.size()*0.60);
      for (std::uint32_t tr = 0; tr < good_threshold; ++tr) {
        auto& good_triplet = tripletsWithPlanes[tr];
        for (std::uint32_t tr2 = std::max(tr+1,bad_threshold); tr2 < tripletsWithPlanes.size(); ++tr2) {
          auto& bad_triplet = tripletsWithPlanes[tr2];
          if((good_triplet.nearSP==bad_triplet.nearSP && 
              (good_triplet.middleSP==bad_triplet.middleSP || good_triplet.farSP==bad_triplet.farSP)) ||
              (good_triplet.middleSP==bad_triplet.middleSP && good_triplet.farSP==bad_triplet.farSP)) {
            //ACTS_INFO("good triplet tr "<<tr<<" has overlap with bad triplet tr2 "<<tr2);
            const Acts::Vector3& abg = tripletsWithPlanes[tr2].getPlaneABG();
            const Acts::ActsScalar& delta = tripletsWithPlanes[tr2].getPlaneDelta();
            // remove this triplet from A and B
            A -= 2. * (abg * abg.transpose());
            B += 2. * delta * abg;
            tripletsWithPlanes.erase(tripletsWithPlanes.begin()+tr2,tripletsWithPlanes.begin()+tr2+1);
            --tr2;
            ++remove_counter;
          }
        }
      }
      ACTS_INFO("removed "<<remove_counter<<" triplets with good-bad triplets / 2");
*/
/*
    }
  }
  ACTS_INFO("B-size of Acts::Vector3 "<<sizeof(Acts::Vector3)<<", "<<sizeof(std::pair<Acts::Vector3, Acts::ActsScalar>)<<"; tripletsWithRays pairs "<<sizeof(std::pair<std::pair<Acts::Vector3, Acts::ActsScalar>, Acts::ActsScalar>)<<", size "<<tripletsWithPlanes.size());
  return vtx;
}
template <typename spacepoint_t>
Acts::Vector3
Acts::SingleSeedVertexFinder<spacepoint_t>::findClosestPointFromRays(
    std::vector<typename Acts::SingleSeedVertexFinder<spacepoint_t>::Triplet>&
        tripletsWithRays) const {
  // 1. define function f = sum over all triplets [distance from an unknown
  // point
  //    (x_0,y_0,z_0) to the ray defined by the triplet]
  // 2. find minimum of "f" by partial derivations over x_0, y_0, and z_0
  // 3. each derivation has parts linearly depending on x_0, y_0, and z_0
  //    (will fill A[][3]) or to nothing (will fill B[])
  // 4. solve A*(x_0,y_0,z_0) = B
  Acts::Vector3 vtx = Acts::Vector3::Zero();
  Acts::Vector3 vtxPrev{m_cfg.rMaxFar, m_cfg.rMaxFar, m_cfg.maxAbsZ};
  // elements of the linear equations to solve
  Acts::SquareMatrix3 A = Acts::SquareMatrix3::Identity() * 2. * tripletsWithRays.size();
  Acts::Vector3 B = Acts::Vector3::Zero();
  for (const auto& triplet : tripletsWithRays) {
    // use ray saved from earlier
    const Acts::Vector3& startPoint=triplet.getStartPoint();
    const Acts::Vector3& direction=triplet.getDirection();
    A -= 2. * (direction * direction.transpose());
    B += -2. * direction * (direction.dot(startPoint)) + 2. * startPoint;
  }
  for (std::uint32_t iter = 0; iter <= m_cfg.maxIterations; iter++) {
    // new vertex position
    vtx = A.lu().solve(B);
    Acts::Vector3 vtxDiff = vtx - vtxPrev;
    ACTS_INFO(iter<<": RAY vtx = "<<vtx[0]<<", "<<vtx[1]<<", "<<vtx[2]<<"; vtxDiff = "<<vtxDiff[0]<<", "<<vtxDiff[1]<<", "<<vtxDiff[2]<<"; there are "<<tripletsWithRays.size()<<" triplets");
    if (vtxDiff.norm() < m_cfg.minVtxShift) {
      // difference between the new vertex and the old vertex is not so large
      break;
    }
    if (iter != m_cfg.maxIterations) {
      // is not the last iteration
      vtxPrev = vtx;
      for (auto& triplet : tripletsWithRays) {
        const Acts::Vector3& startPoint=triplet.getStartPoint();
        const Acts::Vector3& direction=triplet.getDirection();
        const Acts::ActsScalar distance = (vtx - startPoint).cross(direction).norm();
        triplet.setDistance(distance);
      }
      std::sort(tripletsWithRays.begin(), tripletsWithRays.end(),
                [](const auto& lhs, const auto& rhs) {
                  return lhs.getDistance() < rhs.getDistance();
                });
      std::uint32_t threshold = (std::uint32_t)(tripletsWithRays.size() *
                                                (1. - m_cfg.removeFraction));
      for (std::uint32_t tr = threshold + 1; tr < tripletsWithRays.size();
           ++tr) {
        const Acts::Vector3& startPoint = tripletsWithRays[tr].getStartPoint();
        const Acts::Vector3& direction  = tripletsWithRays[tr].getDirection();
        // remove this triplet from A and B
        A -= Acts::SquareMatrix3::Identity() * 2.;
        A += 2. * (direction * direction.transpose());
        B -= -2. * direction * (direction.dot(startPoint)) + 2. * startPoint;
      }
      // remove all excessive triplets
      if(threshold + 1 < tripletsWithRays.size()) {
        ACTS_INFO("removed "<<tripletsWithRays.size()-threshold-1<<" triplets with removeFraction");
        tripletsWithRays.erase(tripletsWithRays.begin()+threshold+1,tripletsWithRays.end());
      }
    }
  }
  return vtx;
}
*/