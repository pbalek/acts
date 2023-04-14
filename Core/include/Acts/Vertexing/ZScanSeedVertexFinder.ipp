// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include <cmath>
#include <Eigen/Eigenvalues> 

template <typename spacepoint_t>
Acts::ZScanSeedVertexFinder<spacepoint_t>::ZScanSeedVertexFinder(const Acts::ZScanSeedVertexFinder<spacepoint_t>::Config& cfg, std::unique_ptr<const Logger> lgr)
    : m_cfg(cfg), m_logger(std::move(lgr))
{
    if(std::isnan(cfg.maxXYdeviation))  ACTS_ERROR("value of maxXYdeviation was not initialized");
    if(std::isnan(cfg.maxXYZdeviation)) ACTS_ERROR("value of maxXYZdeviation was not initialized");
    if(std::isnan(cfg.rMinMiddle))      ACTS_ERROR("value of rMinMiddle was not initialized");
    if(std::isnan(cfg.rMaxMiddle))      ACTS_ERROR("value of rMaxMiddle was not initialized");
    if(std::isnan(cfg.numPhiSlices))    ACTS_ERROR("value of numPhiSlices was not initialized");
    if(cfg.numPhiSlices<3) ACTS_INFO("value of numPhiSlices is "<<cfg.numPhiSlices<<", which is less than 3. There will be duplicate triplets.");
    if(std::isnan(cfg.maxZPosition))    ACTS_ERROR("value of maxZPosition was not initialized");
    if(std::isnan(cfg.maxRPosition))    ACTS_ERROR("value of maxRPosition was not initialized");

    // std::transform(cfg.minimalizeWRT.begin(), cfg.minimalizeWRT.end(), cfg.minimalizeWRT.begin(), [](unsigned char c){ return std::tolower(c); });
    if(cfg.minimalizeWRT!="planes" && cfg.minimalizeWRT!="rays")
    {
        ACTS_ERROR("value of minimalizeWRT is "<<cfg.minimalizeWRT<<", allowed values are \"planes\" or \"rays\" ");
    }
}


template <typename spacepoint_t>
Acts::Vector3 Acts::ZScanSeedVertexFinder<spacepoint_t>::findVertex(const std::vector<spacepoint_t>& spacepoints) const
{
    std::vector<std::vector<std::vector<const spacepoint_t*>>> sorted_spacepoints=sortSpacepoints(spacepoints);

    std::vector<Acts::ZScanSeedVertexFinder<spacepoint_t>::Triplet> triplets=findTriplets(sorted_spacepoints);

    // if no valid triplets found
    if(triplets.empty()) return {};

    Acts::Vector3 vtx=Acts::Vector3::Zero();
    if(m_cfg.minimalizeWRT=="planes")
    {
        // find a point closest to all planes defined by the triplets
        vtx=findClosestPointFromPlanes(triplets);
    }
    else if(m_cfg.minimalizeWRT=="rays")
    {
        // find a point closest to all rays fitted through the triplets
        vtx=findClosestPointFromRays(triplets);
    }

    return vtx;
}


template <typename spacepoint_t>
std::vector<std::vector<std::vector<const spacepoint_t*>>> Acts::ZScanSeedVertexFinder<spacepoint_t>::sortSpacepoints(const std::vector<spacepoint_t>& spacepoints) const
{
    std::vector<const spacepoint_t*> helper={};
    std::vector<std::vector<const spacepoint_t*>> near_spacepoints(m_cfg.numPhiSlices,helper), middle_spacepoints(m_cfg.numPhiSlices,helper), far_spacepoints(m_cfg.numPhiSlices,helper);

    for(const auto& sp : spacepoints)
    {
        Acts::ActsScalar phi = detail::radian_pos(std::atan2(sp.y(),sp.x())); 
        int phislice=(int)(phi/(2*M_PI)*m_cfg.numPhiSlices);

        // every input spacepoint is sorted into one and only one subset
        if(sp.r() < m_cfg.rMinMiddle) near_spacepoints.at(phislice).push_back(&sp);
        else if(sp.r() < m_cfg.rMaxMiddle) middle_spacepoints.at(phislice).push_back(&sp);
        else far_spacepoints.at(phislice).push_back(&sp);
    }

    std::vector<std::vector<std::vector<const spacepoint_t*>>> sorted_spacepoints={near_spacepoints,middle_spacepoints,far_spacepoints};

    return sorted_spacepoints;
}


template <typename spacepoint_t>
std::vector<typename Acts::ZScanSeedVertexFinder<spacepoint_t>::Triplet> Acts::ZScanSeedVertexFinder<spacepoint_t>::findTriplets(const std::vector<std::vector<std::vector<const spacepoint_t*>>>& sorted_spacepoints) const
{
    std::vector<Acts::ZScanSeedVertexFinder<spacepoint_t>::Triplet> triplets;

    for(int nearphi=0;nearphi<m_cfg.numPhiSlices;++nearphi)
    {
        for(const auto& near_sp : sorted_spacepoints.at(0).at(nearphi))
        {
            Acts::ActsScalar phiA = std::atan2(near_sp->y(),near_sp->x());
            for(int middlephi_h=nearphi-1;middlephi_h<=nearphi+1;++middlephi_h)
            {
                int middlephi=(middlephi_h+m_cfg.numPhiSlices)%m_cfg.numPhiSlices;
                for(const auto& middle_sp : sorted_spacepoints.at(1).at(middlephi))
                {
                    Acts::ActsScalar phiB = std::atan2(middle_sp->y(),middle_sp->x());
                    Acts::ActsScalar delta_phiAB=detail::difference_periodic(phiA,phiB,2*M_PI);
                    if(std::abs(delta_phiAB) > m_cfg.maxPhideviation) continue;

                    for(int farphi_h=middlephi-1;farphi_h<=middlephi+1;++farphi_h)
                    {
                        int farphi=(farphi_h+m_cfg.numPhiSlices)%m_cfg.numPhiSlices;
                        for(const auto& far_sp : sorted_spacepoints.at(2).at(farphi))
                        {
                            Acts::ActsScalar phiC = std::atan2(far_sp->y(),far_sp->x());
                            Acts::ActsScalar delta_phiBC=detail::difference_periodic(phiB,phiC,2*M_PI);
                            if(std::abs(delta_phiBC) > m_cfg.maxPhideviation) continue;
                            
                            Acts::ZScanSeedVertexFinder<spacepoint_t>::Triplet tr(near_sp, middle_sp, far_sp);
                            if(isTripletValid(tr)) triplets.push_back(tr);
                        }
                    }
                }
            }
        }
    }

    return triplets;
}


template <typename spacepoint_t>
bool Acts::ZScanSeedVertexFinder<spacepoint_t>::isTripletValid(Acts::ZScanSeedVertexFinder<spacepoint_t>::Triplet& triplet) const
{    
    // slope for near+middle spacepoints
    Acts::ActsScalar alpha1 = std::atan2(triplet.a->y()-triplet.b->y(),triplet.a->x()-triplet.b->x()); 
    // slope for middle+far spacepoints
    Acts::ActsScalar alpha2 = std::atan2(triplet.b->y()-triplet.c->y(),triplet.b->x()-triplet.c->x()); 
    // these two slopes shouldn't be too different
    Acts::ActsScalar delta_alpha=detail::difference_periodic(alpha1,alpha2,2*M_PI);
    if(std::abs(delta_alpha) > m_cfg.maxXYdeviation) return false;

    // near-middle ray
    Acts::Vector3 ab{triplet.a->x()-triplet.b->x(), triplet.a->y()-triplet.b->y(), triplet.a->z()-triplet.b->z()};
    // middle-far ray
    Acts::Vector3 bc{triplet.b->x()-triplet.c->x(), triplet.b->y()-triplet.c->y(), triplet.b->z()-triplet.c->z()};
    // dot product of these two
    Acts::ActsScalar costheta=(ab.dot(bc))/(std::sqrt(ab.dot(ab))*std::sqrt(bc.dot(bc)));
    Acts::ActsScalar theta = std::acos(costheta);
    if(theta>m_cfg.maxXYZdeviation) return false;

    // reject the ray if it doesn't come close to the z-axis
    Acts::Ray3D ray = makeRayFromTriplet(triplet);
    const Acts::Vector3& start_point = ray.origin();
    const Acts::Vector3& direction = ray.dir();
    // norm to z-axis and to the ray
    Acts::Vector3 norm{-1.*direction[1], 1.*direction[0], 0};
    // nearest distance from the ray to z-axis
    Acts::ActsScalar dist = std::fabs(start_point.dot(norm))/std::sqrt(norm.dot(norm));
    if(dist>m_cfg.maxRPosition) return false;

    // cross product of direction and norm
    Acts::Vector3 direction_x_n = direction.cross(norm);
    // z coordinate of the nearest distance from the ray to z-axis
    Acts::ActsScalar zdist = direction_x_n.dot(start_point)/(norm.dot(norm));
    if(std::fabs(zdist)>m_cfg.maxZPosition) return false;
    
    if(m_cfg.minimalizeWRT=="rays")
    {
        // save for later
        triplet.ray=ray;
    }

    return true;
}


template <typename spacepoint_t>
std::pair<Acts::Vector3,Acts::ActsScalar> Acts::ZScanSeedVertexFinder<spacepoint_t>::makePlaneFromTriplet(const Acts::ZScanSeedVertexFinder<spacepoint_t>::Triplet triplet) const
{
    Acts::Vector3 a{triplet.a->x(), triplet.a->y(), triplet.a->z()};
    Acts::Vector3 b{triplet.b->x(), triplet.b->y(), triplet.b->z()};
    Acts::Vector3 c{triplet.c->x(), triplet.c->y(), triplet.c->z()};

    Acts::Vector3 ba=b-a, ca=c-a;

    Acts::Vector3 abg = ba.cross(ca);
    Acts::ActsScalar delta = -1.*abg.dot(a);

    return {abg,delta};
}


template <typename spacepoint_t>
Acts::Vector3 Acts::ZScanSeedVertexFinder<spacepoint_t>::findClosestPointFromPlanes(const std::vector<Acts::ZScanSeedVertexFinder<spacepoint_t>::Triplet>& triplets) const
{
    // define function f = sum over all triplets [distance from an unknown point (x_0,y_0,z_0) to the plane defined by the triplet]
    // find minimum of "f" by partial derivations over x_0, y_0, and z_0
    // each derivation has parts lineary depending on x_0, y_0, and z_0 (will fill A[deriv][3]) or to nothing (will fill B[deriv])
    // solve A*(x_0,y_0,z_0) = B
    
    // elements of the linear equations to solve
    Acts::SymMatrix3 A=Acts::SymMatrix3::Zero();
    Acts::Vector3 B=Acts::Vector3::Zero();

    for(const auto& triplet : triplets)
    {
        auto [abg,delta] = makePlaneFromTriplet(triplet);
        Acts::ActsScalar norm=1./(abg.dot(abg));

        A+=2.*norm*(abg*abg.transpose());
        B-=2.*norm*delta*abg;
    }

    Acts::Vector3 vtx = A.lu().solve(B);

    return vtx;
}


template <typename spacepoint_t>
Acts::Ray3D Acts::ZScanSeedVertexFinder<spacepoint_t>::makeRayFromTriplet(const Acts::ZScanSeedVertexFinder<spacepoint_t>::Triplet triplet) const
{
    Acts::SymMatrix3 mat;
    mat.row(0) = Acts::Vector3(triplet.a->x(), triplet.a->y(), triplet.a->z());
    mat.row(1) = Acts::Vector3(triplet.b->x(), triplet.b->y(), triplet.b->z());
    mat.row(2) = Acts::Vector3(triplet.c->x(), triplet.c->y(), triplet.c->z());

    Acts::Vector3 mean = mat.colwise().mean();
    Acts::SymMatrix3 cov = (mat.rowwise()-mean.transpose()).transpose() * (mat.rowwise()-mean.transpose()) / 3.;

    // "cov" is self-adjoint matrix
    Eigen::SelfAdjointEigenSolver<Acts::SymMatrix3> saes(cov);
    // eigenvalues are sorted in increasing order
    Acts::Vector3 eivec = saes.eigenvectors().col(2);

    return {mean,eivec};
}


template <typename spacepoint_t>
Acts::Vector3 Acts::ZScanSeedVertexFinder<spacepoint_t>::findClosestPointFromRays(const std::vector<Acts::ZScanSeedVertexFinder<spacepoint_t>::Triplet>& triplets) const
{
    // define function f = sum over all triplets [distance from an unknown point (x_0,y_0,z_0) to the ray defined by the triplet]
    // find minimum of "f" by partial derivations over x_0, y_0, and z_0
    // each derivation has parts lineary depending on x_0, y_0, and z_0 (will fill A[][3]) or to nothing (will fill B[])
    // solve A*(x_0,y_0,z_0) = B
    
    // elements of the linear equations to solve
    Acts::SymMatrix3 A=Acts::SymMatrix3::Identity()*2.*triplets.size();
    Acts::Vector3 B=Acts::Vector3::Zero();

    for(const auto& triplet : triplets)
    {
        // use ray saved from earlier
        const Acts::Vector3& start_point = triplet.ray.origin();
        const Acts::Vector3& direction = triplet.ray.dir();

        A-=2.*(direction*direction.transpose());
        B+=-2.*direction*(direction.dot(start_point))+2.*start_point;
    }

    Acts::Vector3 vtx = A.lu().solve(B);

    return vtx;
}
