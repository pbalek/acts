// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include <cmath>
#include <utility>
#include <iostream>

template <typename spacepoint_t>
Acts::ZScanSeedVertexFinder<spacepoint_t>::ZScanSeedVertexFinder(const Acts::ZScanSeedVertexFinder<spacepoint_t>::Config& cfg, std::unique_ptr<const Logger> lgr)
    : m_cfg(cfg), m_logger(std::move(lgr))
{
    if(std::isnan(cfg.maxZRdeviation)) ACTS_ERROR("value of maxZRdeviation was not initialized");
    if(std::isnan(cfg.maxXYdeviation)) ACTS_ERROR("value of maxXYdeviation was not initialized");
    if(std::isnan(cfg.rMinMiddle))     ACTS_ERROR("value of rMinMiddle was not initialized");
    if(std::isnan(cfg.rMaxMiddle))     ACTS_ERROR("value of rMaxMiddle was not initialized");
    if(std::isnan(cfg.zBinSize))       ACTS_ERROR("value of zBinSize was not initialized");
    if(std::isnan(cfg.maxZPosition))   ACTS_ERROR("value of maxZPosition was not initialized");
}


template <typename spacepoint_t>
std::vector<float> Acts::ZScanSeedVertexFinder<spacepoint_t>::findVertex(const std::vector<spacepoint_t>& spacepoints) const
{
    // sort all spacepoints to new vectors
    std::vector<std::vector<const spacepoint_t*>> sorted_spacepoints=sortSpacepoints(spacepoints);

    // find all valid triplets
    std::vector<Acts::ZScanSeedVertexFinder<spacepoint_t>::Triplet> triplets=findTriplets(sorted_spacepoints);

    // if no valid triplets found
    if(triplets.empty()) return {};

    // create a vector pretending to be a histogram of the estimated origins of the 
    std::vector<int> hist=makeZHist(triplets);

    // find a peak in the provided histogram; may get confused if the peak is not clear
    float peak=findZPeak(hist);

    // position along z-axis; takes "zBinSize" into account
    return {peak};
}


template <typename spacepoint_t>
std::vector<std::vector<const spacepoint_t*>> Acts::ZScanSeedVertexFinder<spacepoint_t>::sortSpacepoints(const std::vector<spacepoint_t>& spacepoints) const
{
    std::vector<const spacepoint_t*> near_spacepoints, middle_spacepoints, far_spacepoints;

    for(const auto& sp : spacepoints)
    {
        // every input spacepoint is sorted into one and only one subset
        if(sp.r() < m_cfg.rMinMiddle) near_spacepoints.push_back(&sp);
        else if(sp.r() < m_cfg.rMaxMiddle) middle_spacepoints.push_back(&sp);
        else far_spacepoints.push_back(&sp);
    }

    std::vector<std::vector<const spacepoint_t*>> sorted_spacepoints={near_spacepoints,middle_spacepoints,far_spacepoints};

    return sorted_spacepoints;
}


template <typename spacepoint_t>
std::vector<typename Acts::ZScanSeedVertexFinder<spacepoint_t>::Triplet> Acts::ZScanSeedVertexFinder<spacepoint_t>::findTriplets(const std::vector<std::vector<const spacepoint_t*>>& sorted_spacepoints) const
{
    std::vector<Acts::ZScanSeedVertexFinder<spacepoint_t>::Triplet> triplets;

    for(const auto& near_sp : sorted_spacepoints.at(0))
    {
        for(const auto& middle_sp : sorted_spacepoints.at(1))
        {
            for(const auto& far_sp : sorted_spacepoints.at(2))
            {
                Acts::ZScanSeedVertexFinder<spacepoint_t>::Triplet tr(near_sp, middle_sp, far_sp);
                if(isTripletValid(tr)) triplets.push_back(tr);
            }
        }
    }

    return triplets;
}


template <typename spacepoint_t>
bool Acts::ZScanSeedVertexFinder<spacepoint_t>::isTripletValid(const Acts::ZScanSeedVertexFinder<spacepoint_t>::Triplet triplet) const
{
    // TODO: check all quadrants
   
    float a1 = (triplet.a->y() - triplet.b->y())/(triplet.a->x() - triplet.b->x());
    float alpha1 = std::atan(a1); // slope for near+middle spacepoints

    float a2 = (triplet.b->y() - triplet.c->y())/(triplet.b->x() - triplet.c->x());
    float alpha2 = std::atan(a2); // slope for middle+far spacepoints

    // these two slopes shouldn't be too different
    if(std::abs(alpha2-alpha1) > m_cfg.maxXYdeviation) return false;

    float b1 = (triplet.a->r() - triplet.b->r())/(triplet.a->z() - triplet.b->z());
    float beta1 = std::atan(b1); // slope for near+middle spacepoints

    float b2 = (triplet.b->r() - triplet.c->r())/(triplet.b->z() - triplet.c->z());
    float beta2 = std::atan(b2); // slope for middle+far spacepoints

    // these two slopes shouldn't be too different
    if(std::abs(beta2-beta1) > m_cfg.maxZRdeviation) return false;
    
    return true;
}


template <typename spacepoint_t>
std::vector<int> Acts::ZScanSeedVertexFinder<spacepoint_t>::makeZHist(const std::vector<Acts::ZScanSeedVertexFinder<spacepoint_t>::Triplet>& triplets) const
{
    std::vector<int> hist;
    // convention for the bin numbers:
    //   -------------|-------------
    //                0          -> z
    // ..|11|9|7|5|3|1|0|2|4|6|8|10|12|..
    hist.resize(200,0);

    for(const auto& triplet : triplets)
    {
        // TODO: check if this is valid
        
        // fit the 3 spacepoints with a straight line: r = slope*z + cons
        float Sz=0.,Sr=0.,Szz=0.,Szr=0.;
        std::vector<const spacepoint_t*> spacepoints = {triplet.a,triplet.b,triplet.c};
        for(const auto& sp : spacepoints)
        {
            Sz+=sp->z();
            Sr+=sp->r();
            Szz+=sp->z()*sp->z();
            Szr+=sp->z()*sp->r();
        }
        float slope=(3*Szr - Sz*Sr)/(3*Szz - Sz*Sz);
        float cons=(Sr-slope*Sz)/3.;

        // estimated position of the vertex for this triplet
        float z=-1.*cons/slope;

        // ignore possible vertices that are too far; this also prevents the histogram to grow too much
        if(std::fabs(z)>m_cfg.maxZPosition) continue;

        // std::cout<<" spacepoint a "<<triplet.a->r()<<" "<<triplet.a->z()<<std::endl;
        // std::cout<<" spacepoint b "<<triplet.b->r()<<" "<<triplet.b->z()<<std::endl;
        // std::cout<<" spacepoint c "<<triplet.c->r()<<" "<<triplet.c->z()<<std::endl;
        ACTS_DEBUG("straight line fit is r = "<<slope<<"*z + "<<cons);

        // desired histogram bin for this triplet
        unsigned int zbin=2*(unsigned int)(std::fabs(z)/m_cfg.zBinSize);
        if(z<0) ++zbin;
        
        std::cout<<"this triplet has z = "<<z<<" and will go to the bin "<<zbin<<std::endl;

        // extend the histogram with zeros, if needed 
        while(zbin > hist.size()-1)
        {
            hist.resize(2*hist.size(),0);
            ACTS_DEBUG("new hist size = "<<hist.size());
        }

        ++hist.at(zbin);
    }

    return hist;
}


template <typename spacepoint_t>
float Acts::ZScanSeedVertexFinder<spacepoint_t>::findZPeak(const std::vector<int>& hist) const
{
    // find maximum
    int maxh=0;
    for(unsigned int h=0;h<hist.size();h++)
    {
        if(hist.at(h) > hist.at(maxh)) maxh=h;
    }

    // look around the maximum for a better estimation of the peak
    std::vector<int> hpeak{maxh-2,maxh,maxh+2};
    // fix for bins around zero
    if(maxh==0) hpeak.at(0)=1;
    if(maxh==1) hpeak.at(0)=0;

    float zsum=0.,zpos=0.;
    for(auto h : hpeak)
    {
        // position of the edge of the bin closer to zero
        float z=m_cfg.zBinSize*(h/2);
        // increase/decrease the position by a half bin size, so it corresponds to the middle of the bin
        if(h%2) z-=m_cfg.zBinSize/2.;
        else    z+=m_cfg.zBinSize/2.;

        // make weighted average
        zpos+=z*hist.at(h);
        zsum+=hist.at(h);
    }

    return zpos/zsum;
}
