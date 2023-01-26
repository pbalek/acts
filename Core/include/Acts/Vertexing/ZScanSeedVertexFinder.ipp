// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename spacepoint_t>
Acts::ZScanSeedVertexFinder<spacepoint_t>::ZScanSeedVertexFinder(Acts::ZScanSeedVertexFinder<spacepoint_t>::Config& cfg)
    : m_cfg(cfg)
{
    if(std::isnan(cfg.maxZRdeviation))
    {
        throw std::runtime_error("value of maxZRdeviation was not initialized");
    }
    if(std::isnan(cfg.maxXYdeviation))
    {
        throw std::runtime_error("value of maxXYdeviation was not initialized");
    }
}


template <typename spacepoint_t>
std::vector<float> Acts::ZScanSeedVertexFinder<spacepoint_t>::findVertex(const std::vector<spacepoint_t>& spacepoints) const
{
    std::vector<std::vector<spacepoint_t>> sorted_spacepoints=sortSpacepoints(spacepoints);

    std::vector<Triplet> triplets=getTriplets(sorted_spacepoints);

    std::vector<int> hist=makeZHist(triplets);

    float peak=findZPeak(hist);

    return {peak};
}


template <typename spacepoint_t>
std::vector<std::vector<spacepoint_t>> Acts::ZScanSeedVertexFinder<spacepoint_t>::sortSpacepoints(const std::vector<spacepoint_t>& spacepoints) const
{
    std::vector<spacepoint_t> near_spacepoints, middle_spacepoints, far_spacepoints;

    for(auto sp : spacepoints)
    {
        if(sp.radius() < m_cfg.rMinMiddle) near_spacepoints.push_back(sp);
        else if(sp.radius() < m_cfg.rMaxMiddle) middle_spacepoints.push_back(sp);
        else far_spacepoints.push_back(sp);
    }

    std::vector<std::vector<spacepoint_t>> sorted_spacepoints={near_spacepoints,middle_spacepoints,far_spacepoints};

    return sorted_spacepoints;
}


template <typename spacepoint_t>
std::vector<Triplet> Acts::ZScanSeedVertexFinder<spacepoint_t>::findTriplets(const std::vector<std::vector<spacepoint_t>>& sorted_spacepoints) const
{
    std::vector<Triplet> triplets;

    for(auto near_sp : sorted_spacepoint.at(0))
    {
        for(auto middle_sp : sorted_spacepoint.at(1))
        {
            for(auto far_sp : sorted_spacepoint.at(2))
            {
                Triplet tr(near_sp, middle_sp, far_sp);
                if(isValidTriplet(tr)) triplets.push_back(tr);
            }
        }
    }

    return triplets;
}


template <typename spacepoint_t>
bool Acts::ZScanSeedVertexFinder<spacepoint_t>::isTripletValid(const Triplet triplet) const
{
    // TODO: check all quadrants
   
    float a1 = (triplet.a.y() - triplet.b.y())/(triplet.a.x() - triplet.b.x());
    float alpha1 = std::atan(a1);

    float a2 = (triplet.b.y() - triplet.c.y())/(triplet.b.x() - triplet.c.x())
    float alpha2 = std::atan(a2);

    if(std::abs(alpha2-alpha1) > m_cfg.maxXYdeviation) return false;

    float b1 = (triplet.a.r() - triplet.b.r())/(triplet.a.z() - triplet.b.z());
    float beta1 = std::atan(b1);

    float b2 = (triplet.b.r() - triplet.c.r())/(triplet.b.z() - triplet.c.z())
    float beta2 = std::atan(b2);

    if(std::abs(beta2-beta1) > m_cfg.maxZRdeviation) return false;
    
    return true;
}


template <typename spacepoint_t>
std::vector<int> Acts::ZScanSeedVertexFinder<spacepoint_t>::makeZHist(const std::vector<Triplet>& triplets) const
{
    std::vector<int> hist;
    //   -------------|-------------
    //                0          -> z
    // ..|11|9|7|5|3|1|0|2|4|6|8|10|12|..
    hist.reserve(200,0);

    for(auto triplet : triplets)
    {
        // TODO: check if this is valid
        
        float Sz=0.,Sr=0.,Szz=0.,Szr=0.;
        Sz=a.z();
        Sr=a.radius();
        Szz=a.z()*a.z();
        Szr=a.z()*a.radius();

        Sz+=b.z();
        Sr+=b.radius();
        Szz+=b.z()*b.z();
        Szr+=b.z()*b.radius();

        Sz+=c.z();
        Sr+=c.radius();
        Szz+=c.z()*c.z();
        Szr+=c.z()*c.radius();

        float slope=(3*Szr - Sz*Sr)/(3*Szz - Sz*Sz);
        float cons=(Sr-slope*Sz)/3.;

        float z=-1.*slope/cons;

        int zbin=2*int(std::abs(z)/m_cgf.zBinSize)
        if(z<0) zbin-=1;

        if(zbin > hist.size()-1) hist.resize(2*hist.size(),0);

        ++hist.at(zbin);
    }

    return hist;
}


template <typename spacepoint_t>
float Acts::ZScanSeedVertexFinder<spacepoint_t>::findZPeak(std::vector<float>& hist) const;
{
    int maxh=0;

    for(unsigned int h=0;h<hist.size();h++)
    {
        if(hist.at(h) > hist.at(maxh)) maxh=h;
    }

    float z=m_cgf.zBinSize*(maxh/2);
    if(z%2) z+=m_cgf.zBinSize/2.;
    else    z-=m_cgf.zBinSize/2.;

    return z;
}
