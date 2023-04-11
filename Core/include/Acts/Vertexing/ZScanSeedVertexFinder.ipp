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

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 

template <typename spacepoint_t>
Acts::ZScanSeedVertexFinder<spacepoint_t>::ZScanSeedVertexFinder(const Acts::ZScanSeedVertexFinder<spacepoint_t>::Config& cfg, std::unique_ptr<const Logger> lgr)
    : m_cfg(cfg), m_logger(std::move(lgr))
{
    if(std::isnan(cfg.maxZRdeviation)) ACTS_ERROR("value of maxZRdeviation was not initialized");
    if(std::isnan(cfg.maxXYdeviation)) ACTS_ERROR("value of maxXYdeviation was not initialized");
    if(std::isnan(cfg.rMinMiddle))     ACTS_ERROR("value of rMinMiddle was not initialized");
    if(std::isnan(cfg.rMaxMiddle))     ACTS_ERROR("value of rMaxMiddle was not initialized");
    if(std::isnan(cfg.minZRslope))     ACTS_ERROR("value of minZRslope was not initialized");
    if(std::isnan(cfg.zBinSize))       ACTS_ERROR("value of zBinSize was not initialized");
    if(std::isnan(cfg.maxZPosition))   ACTS_ERROR("value of maxZPosition was not initialized");
}


template <typename spacepoint_t>
std::vector<std::vector<float>> Acts::ZScanSeedVertexFinder<spacepoint_t>::findVertex(const std::vector<spacepoint_t>& spacepoints) const
{
    // sort all spacepoints to new vectors
    std::vector<std::vector<const spacepoint_t*>> sorted_spacepoints=sortSpacepoints(spacepoints);

    // find all valid triplets
    std::vector<Acts::ZScanSeedVertexFinder<spacepoint_t>::Triplet> triplets=findTriplets(sorted_spacepoints);

    // if no valid triplets found
    if(triplets.empty()) return {};

    // find a point closest to all planes defined by the triplets
    std::vector<float> vtx=findClosestPointFromPlanes(triplets);

    // find a point closest to all lines fitted through the triplets
    std::vector<float> vtx2=findClosestPointFromLines(triplets);

    std::vector<std::vector<float>> vtx3=findClosestPointFromZRLines(triplets);

    // create a vector pretending to be a histogram of the estimated origins of the 
    std::vector<int> hist=makeZHist(triplets);
    std::vector<float> hist2={};
    for(auto h : hist) hist2.push_back(1.0*h);

    // find a peak in the provided histogram; may get confused if the peak is not clear
    float peak=findZPeak(hist,0);
    std::vector<float> peakvtx{0.,0.,peak};

    float peak1=findZPeak(hist,1);
    std::vector<float> peakvtx1{0.,0.,peak1};

    float peak2=findZPeak(hist,2);
    std::vector<float> peakvtx2{0.,0.,peak2};

    // position along z-axis; takes "zBinSize" into account
    // return {peakvtx,hist2};

    return {peakvtx,peakvtx1,peakvtx2,vtx,vtx2,vtx3.at(0),vtx3.at(1)};

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

    std::cout<<"have spacepoints: near "<<near_spacepoints.size()<<", middle "<<middle_spacepoints.size()<<", far "<<far_spacepoints.size()<<std::endl;

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
            float phiA = std::atan2(near_sp->y(),near_sp->x());
            float phiB = std::atan2(middle_sp->y(),middle_sp->x());
            float delta_phiAB=std::abs(phiA-phiB);
            if(delta_phiAB > 3.141592654) delta_phiAB = 2*3.141592654-delta_phiAB;
            if(std::abs(delta_phiAB) > m_cfg.maxXYdeviation) return false;
            
            for(const auto& far_sp : sorted_spacepoints.at(2))
            {
                float phiC = std::atan2(far_sp->y(),far_sp->x());
                float delta_phiBC=std::abs(phiB-phiC);
                if(delta_phiBC > 3.141592654) delta_phiBC = 2*3.141592654-delta_phiBC;
                if(std::abs(delta_phiBC) > m_cfg.maxXYdeviation) return false;
                
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
    // // slope for near+middle spacepoints
    // float beta1 = std::atan2(triplet.a->r() - triplet.b->r(),triplet.a->z() - triplet.b->z());
    // if(std::abs(beta1) < m_cfg.minZRslope || (2*3.141592654-std::abs(beta1))<m_cfg.minZRslope) return false;
    // // slope for middle+far spacepoints
    // float beta2 = std::atan2(triplet.b->r() - triplet.c->r(),triplet.b->z() - triplet.c->z()); 
    // // these two slopes shouldn't be too different
    // float delta_beta=std::abs(beta2-beta1);
    // if(delta_beta > 3.141592654) delta_beta = 2*3.141592654-delta_beta;
    // if(std::abs(delta_beta) > m_cfg.maxZRdeviation) return false;

    // slope for near+middle spacepoints
    float alpha1 = std::atan2(triplet.a->y() - triplet.b->y(),triplet.a->x() - triplet.b->x()); 
    // slope for middle+far spacepoints
    float alpha2 = std::atan2(triplet.b->y() - triplet.c->y(),triplet.b->x() - triplet.c->x()); 
    // these two slopes shouldn't be too different
    float delta_alpha=std::abs(alpha2-alpha1);
    if(delta_alpha > 3.141592654) delta_alpha = 2*3.141592654-delta_alpha;
    if(std::abs(delta_alpha) > m_cfg.maxXYdeviation) return false;

    float ab_x=(triplet.a->x() - triplet.b->x());
    float ab_y=(triplet.a->y() - triplet.b->y());
    float ab_z=(triplet.a->z() - triplet.b->z());

    float bc_x=(triplet.b->x() - triplet.c->x());
    float bc_y=(triplet.b->y() - triplet.c->y());
    float bc_z=(triplet.b->z() - triplet.c->z());

    float costheta=(ab_x*bc_x + ab_y*bc_y + ab_z*bc_z)/(std::sqrt(ab_x*ab_x + ab_y*ab_y + ab_z*ab_z)*std::sqrt(bc_x*bc_x + bc_y*bc_y + bc_z*bc_z));
    float theta = std::acos(costheta);
    if(theta>0.1) return false;


    // reject the line if it doesn't come close to the z-axis
    std::array<double,6> line = makeLineFromTriplet(triplet);
    std::vector<double> start_point{line[0],line[1],line[2]};
    std::vector<double> direction{line[3],line[4],line[5]};

    std::vector<double> n{-1.*direction[1], 1.*direction[0], 0};
    double d = std::fabs( start_point[0]*n[0] + start_point[1]*n[1])/std::sqrt(n[0]*n[0]+n[1]*n[1]);
    if(d>10.) return false;

    std::vector<double> direction_x_n{direction[1]*n[2]-direction[2]*n[1], direction[2]*n[0]-direction[0]*n[1], direction[0]*n[1]-direction[1]*n[0]};
    double t2 = (direction_x_n[0]*start_point[0]+direction_x_n[1]*start_point[1]+direction_x_n[2]*start_point[2])/(n[0]*n[0]+n[1]*n[1]);

    if(std::fabs(t2)>m_cfg.maxZPosition) return false;

    // if(theta<std::abs(delta_alpha) || theta<std::abs(delta_beta))
    // {
    //     std::cout<<"theta = "<<theta<<"; delta xy = "<<delta_alpha<<", delta rz = "<<delta_beta<<std::endl;
    //     std::cout<<"  a = "<<triplet.a->x()<<", "<<triplet.a->y()<<", "<<triplet.a->z()<<std::endl;
    //     std::cout<<"  b = "<<triplet.b->x()<<", "<<triplet.b->y()<<", "<<triplet.b->z()<<std::endl;
    //     std::cout<<"  c = "<<triplet.c->x()<<", "<<triplet.c->y()<<", "<<triplet.c->z()<<std::endl;
    // }

    return true;
}


template <typename spacepoint_t>
std::array<double,4> Acts::ZScanSeedVertexFinder<spacepoint_t>::makePlaneFromTriplet(const Acts::ZScanSeedVertexFinder<spacepoint_t>::Triplet triplet) const
{
    double ax=triplet.a->x(), ay=triplet.a->y(), az=triplet.a->z();
    double bx=triplet.b->x(), by=triplet.b->y(), bz=triplet.b->z();
    double cx=triplet.c->x(), cy=triplet.c->y(), cz=triplet.c->z();

    double alpha = (by-ay)*(cz-az) - (cy-ay)*(bz-az);
    double beta  = (bz-az)*(cx-ax) - (cz-az)*(bx-ax);
    double gamma = (bx-ax)*(cy-ay) - (cx-ax)*(by-ay);
    double delta = -alpha*ax - beta*ay - gamma*az; 

    return {alpha,beta,gamma,delta};
}


template <typename spacepoint_t>
std::vector<float> Acts::ZScanSeedVertexFinder<spacepoint_t>::findClosestPointFromPlanes(const std::vector<Acts::ZScanSeedVertexFinder<spacepoint_t>::Triplet>& triplets) const
{
    // define function f = sum over all triplets [distance from an unknown point (x_0,y_0,z_0) to the plane defined by the triplet]
    // find minimum of "f" by partial derivations over x_0, y_0, and z_0
    // each derivation has parts lineary depending on x_0, y_0, and z_0 (will fill A[][3]) or to nothing (will fill B[])
    // solve A*(x_0,y_0,z_0) = B
    
    // elements of the linear equations to solve
    Eigen::Matrix3d A=Eigen::Matrix3d::Zero(3,3);
    Eigen::Vector3d B=Eigen::Vector3d::Zero(3);

    // int counter=0;
    for(const auto& triplet : triplets)
    {
        // auto [alpha,beta,gamma,delta] = makePlaneFromTriplet(triplet);
        std::array<double,4> plane = makePlaneFromTriplet(triplet);
        double norm=1./(plane[0]*plane[0]+plane[1]*plane[1]+plane[2]*plane[2]);

        for(int deriv=0;deriv<3;++deriv)
        {
            for(int depend=0; depend<3;++depend)
            {
                // part of df/d{deriv} that depends on {depend}
                A(deriv,depend)+=norm*(2*plane[deriv]*plane[depend]);
            }

            B(deriv)+=norm*(-2*plane[deriv]*plane[3]);
        }

        // ++counter;
        // if(counter==3) break;
    }

    Eigen::Vector3d x = A.lu().solve(B);

    return {(float)x(0),(float)x(1),(float)x(2)};
}


template <typename spacepoint_t>
std::array<double,6> Acts::ZScanSeedVertexFinder<spacepoint_t>::makeLineFromTriplet(const Acts::ZScanSeedVertexFinder<spacepoint_t>::Triplet triplet) const
{
    std::vector<std::vector<double>> triplet_vec{{triplet.a->x(), triplet.a->y(), triplet.a->z()},
                                                 {triplet.b->x(), triplet.b->y(), triplet.b->z()},
                                                 {triplet.c->x(), triplet.c->y(), triplet.c->z()}};

    Eigen::Vector3d mean=Eigen::Vector3d::Zero(3);
    Eigen::Matrix3d corr=Eigen::Matrix3d::Zero(3,3);

    for(int abc=0;abc<3;++abc)
    {
        for(int xyz=0;xyz<3;++xyz)
        {
            mean(xyz)+=triplet_vec.at(abc).at(xyz);
        }

        for(int row=0;row<3;++row)
        {
            for(int col=0;col<3;++col)
            {
                corr(row,col)+=triplet_vec.at(abc).at(row)*triplet_vec.at(abc).at(col);
            }
        }
    }
    mean/=3.;
    corr/=3.;

    Eigen::Matrix3d cov;
    cov << corr(0,0)-mean(0)*mean(0), corr(0,1)-mean(0)*mean(1), corr(0,2)-mean(0)*mean(2),
           corr(0,1)-mean(0)*mean(1), corr(1,1)-mean(1)*mean(1), corr(1,2)-mean(1)*mean(2),
           corr(0,2)-mean(0)*mean(2), corr(1,2)-mean(2)*mean(1), corr(2,2)-mean(2)*mean(2);
    
    Eigen::EigenSolver<Eigen::Matrix3d> es(cov);

    // eigen values are real and non-negative, so "pseudo" is Ok
    Eigen::Matrix3d eivals = es.pseudoEigenvalueMatrix();
    Eigen::Matrix3d eivecs = es.pseudoEigenvectors();

    // find largest eigenvalue and use that eigenvector as the direction
    // eigenvalue matrix is always diagonal
    int maxval;
    if(eivals(1,1)>eivals(0,0)) maxval=1;
    else maxval=0;
    if(eivals(2,2)>eivals(maxval)) maxval=2;

    return {mean(0),mean(1),mean(2),eivecs(0,maxval),eivecs(1,maxval),eivecs(2,maxval)};
}


template <typename spacepoint_t>
std::vector<float> Acts::ZScanSeedVertexFinder<spacepoint_t>::findClosestPointFromLines(const std::vector<Acts::ZScanSeedVertexFinder<spacepoint_t>::Triplet>& triplets) const
{
    // define function f = sum over all triplets [distance from an unknown point (x_0,y_0,z_0) to the line defined by the triplet]
    // find minimum of "f" by partial derivations over x_0, y_0, and z_0
    // each derivation has parts lineary depending on x_0, y_0, and z_0 (will fill A[][3]) or to nothing (will fill B[])
    // solve A*(x_0,y_0,z_0) = B
    
    // elements of the linear equations to solve
    Eigen::Matrix3d A=Eigen::Matrix3d::Zero(3,3);
    Eigen::Vector3d B=Eigen::Vector3d::Zero(3);
    
    for(const auto& triplet : triplets)
    {
        // auto [alpha,beta,gamma,delta] = makePlaneFromTriplet(triplet);
        std::array<double,6> line = makeLineFromTriplet(triplet);
        std::vector<double> start_point{line[0],line[1],line[2]};
        std::vector<double> direction{line[3],line[4],line[5]};

        // reject the line if it doesn't come close to the z-axis
        std::vector<double> n{-1.*direction[1], 1.*direction[0], 0};
        double d = std::fabs( start_point[0]*n[0] + start_point[1]*n[1])/std::sqrt(n[0]*n[0]+n[1]*n[1]);
        if(d>10.) continue;

        for(int deriv=0;deriv<3;++deriv)
        {
            for(int depend=0; depend<3;++depend)
            {
                // part of df/d{deriv} that depends on {depend}
                if(deriv==depend) A(deriv,depend)+=2.*(1.-direction[deriv]*direction[deriv]);
                else              A(deriv,depend)+=-2.*direction[deriv]*direction[depend];

                // part of df/d{deriv} that doesn't depends on anothing
                if(deriv==depend) B(deriv)+=2.*(1.-direction[deriv]*direction[deriv])*start_point[depend];
                else              B(deriv)+=-2*direction[deriv]*direction[depend]*start_point[depend];
            }
        }
    }

    Eigen::Vector3d x = A.lu().solve(B);

    // std::vector<double> wannabevtx{(float)x(0),(float)x(1),(float)x(2)};
    // std::vector<double> truevtx{1.5,0.,20.};

    // float chi2=0.,truevtxchi2=0.;
    // int tripletcounter=0;
    // float maxchi2=0., truemaxchi2=0.;
    // for(const auto& triplet : triplets)
    // {
    //     ++tripletcounter;

    //     // auto [alpha,beta,gamma,delta] = makePlaneFromTriplet(triplet);
    //     std::array<double,6> line = makeLineFromTriplet(triplet);
    //     std::vector<double> start_point{line[0],line[1],line[2]};
    //     std::vector<double> direction{line[3],line[4],line[5]};

    //     // reject the line if it doesn't come close to the z-axis
    //     std::vector<double> n{-1.*direction[1], 1.*direction[0], 0};
    //     double d = std::fabs( start_point[0]*n[0] + start_point[1]*n[1])/std::sqrt(n[0]*n[0]+n[1]*n[1]);

    //     std::vector<double> direction_x_n{direction[1]*n[2]-direction[2]*n[1], direction[2]*n[0]-direction[0]*n[1], direction[0]*n[1]-direction[1]*n[0]};
    //     double t2 = (direction_x_n[0]*start_point[0]+direction_x_n[1]*start_point[1]+direction_x_n[2]*start_point[2])/(n[0]*n[0]+n[1]*n[1]);
    //     if(d>10.) continue;

    //     std::vector<double> vec{wannabevtx[0]-start_point[0],wannabevtx[1]-start_point[1],wannabevtx[2]-start_point[2]};
    //     double dist = pow(vec[1]*direction[2]-vec[2]*direction[1],2) + pow(vec[2]*direction[0]-vec[0]*direction[2],2) + pow(vec[0]*direction[1]-vec[1]*direction[0],2);
    //     if(maxchi2<dist) maxchi2=dist;
    //     chi2+=dist;

    //     std::vector<double> vec2{truevtx[0]-start_point[0],truevtx[1]-start_point[1]};
    //     dist = pow(vec2[1]*direction[2]-vec2[2]*direction[1],2) + pow(vec2[2]*direction[0]-vec2[0]*direction[2],2) + pow(vec2[0]*direction[1]-vec2[1]*direction[0],2);
    //     if(truemaxchi2<dist) truemaxchi2=dist;
    //     truevtxchi2+=dist;

    //     if(dist > 1e4)
    //     {
    //         std::cout<<" . line dealing "<<dist<<" chi2: start "<<line[0]<<","<<line[1]<<","<<line[2]<<", direction "<<line[3]<<","<<line[4]<<","<<line[5]<<"; distance from z = "<<d<<", position at z "<<t2<<std::endl;
    //     }
    // }

    // std::cout<<" .tripletcounter = "<<tripletcounter<<std::endl;
    // std::cout<<" .wannabevtx chi2 = "<<chi2<<", true vtx chi2 = "<<truevtxchi2<<std::endl;
    // std::cout<<" .       max chi2 = "<<maxchi2<<",    max chi2 = "<<truemaxchi2<<std::endl;


    return {(float)x(0),(float)x(1),(float)x(2)};
}

template <typename spacepoint_t>
std::array<double,2> Acts::ZScanSeedVertexFinder<spacepoint_t>::makeZRLineFromTriplet(const Acts::ZScanSeedVertexFinder<spacepoint_t>::Triplet triplet) const
{
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

    return {slope,cons};
}

template <typename spacepoint_t>
std::vector<std::vector<float>> Acts::ZScanSeedVertexFinder<spacepoint_t>::findClosestPointFromZRLines(const std::vector<Acts::ZScanSeedVertexFinder<spacepoint_t>::Triplet>& triplets) const
{
    // define function f = sum over all triplets [distance from an unknown point (x_0,y_0,z_0) to the line defined by the triplet]
    // find minimum of "f" by partial derivations over x_0, y_0, and z_0
    // each derivation has parts lineary depending on x_0, y_0, and z_0 (will fill A[][3]) or to nothing (will fill B[])
    // solve A*(x_0,y_0,z_0) = B
    
    // elements of the linear equations to solve
    Eigen::Matrix2d A=Eigen::Matrix2d::Zero(2,2);
    Eigen::Vector2d B=Eigen::Vector2d::Zero(2);
    
    std::vector<float> zrlines{};

    std::cout<<"have "<<triplets.size()<<" triplets"<<std::endl;

    for(const auto& triplet : triplets)
    {
        // auto [alpha,beta,gamma,delta] = makePlaneFromTriplet(triplet);
        std::array<double,2> line = makeZRLineFromTriplet(triplet);
        // ignore possible vertices that are too far; this also prevents the histogram to grow too much
        if(std::fabs(line[0]) < m_cfg.minZRslope) continue;
        if(std::fabs(line[1]/line[0])>m_cfg.maxZPosition) continue;

        zrlines.push_back(line[0]);
        zrlines.push_back(line[1]);

        std::vector<double> start_point{-line[1]/line[0],0.};
        std::vector<double> direction_norm{-1.*std::sin(std::atan(line[0])),std::cos(std::atan(line[0]))};

        for(int deriv=0;deriv<2;++deriv)
        {
            for(int depend=0; depend<2;++depend)
            {
                // part of df/d{deriv} that depends on {depend}
                A(deriv,depend)+=(2*direction_norm[deriv]*direction_norm[depend]);

                // part of df/d{deriv} that doesn't depends on anothing
                B(deriv)+=(2*direction_norm[deriv]*direction_norm[depend]*start_point[depend]);
            }
        }
    }

    Eigen::Vector2d x = A.lu().solve(B);

    // std::vector<double> wannabevtx{(double)x(0),(double)x(1)};
    // std::vector<double> truevtx{20.,1.5};

    // float chi2=0.,truevtxchi2=0.;
    // int tripletcounter=0;
    // float maxchi2=0., truemaxchi2=0.;
    // for(const auto& triplet : triplets)
    // {
    //     ++tripletcounter;

    //     std::array<double,2> line = makeZRLineFromTriplet(triplet);
    //     // ignore possible vertices that are too far; this also prevents the histogram to grow too much
    //     if(std::fabs(line[0]) < m_cfg.minZRslope) continue;
    //     if(std::fabs(line[1]/line[0])>m_cfg.maxZPosition) continue;

    //     std::vector<double> start_point{-line[1]/line[0], 0.};
    //     std::vector<double> direction_norm{-1.*std::sin(std::atan(line[0])),std::cos(std::atan(line[0]))};
    //     // float weight=std::fabs(line[0]);
    //     // float weight=1.;

    //     std::vector<double> vec{wannabevtx[0]-start_point[0],wannabevtx[1]-start_point[1]};
    //     double dist = vec[0]*direction_norm[0]+vec[1]*direction_norm[1];
    //     if(maxchi2<dist) maxchi2=dist;
    //     chi2+=dist;

    //     std::vector<double> vec2{truevtx[0]-start_point[0],truevtx[1]-start_point[1]};
    //     dist = vec2[0]*direction_norm[0]+vec2[1]*direction_norm[1];
    //     if(truemaxchi2<dist) truemaxchi2=dist;
    //     truevtxchi2+=dist;

    //     if(dist > 1e3)
    //     {
    //         std::cout<<" - line dealing "<<dist<<" chi2: slope "<<line[0]<<", const "<<line[1]<<std::endl;
    //     }
    // }

    // std::cout<<"tripletcounter = "<<tripletcounter<<std::endl;
    // std::cout<<"wannabevtx chi2 = "<<chi2<<", true vtx chi2 = "<<truevtxchi2<<std::endl;
    // std::cout<<"       max chi2 = "<<maxchi2<<",    max chi2 = "<<truemaxchi2<<std::endl;

    return {{(float)x(1),0.,(float)x(0)},zrlines};
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
        // ACTS_DEBUG("straight line fit is r = "<<slope<<"*z + "<<cons);

        // desired histogram bin for this triplet
        unsigned int zbin=2*(unsigned int)(std::fabs(z)/m_cfg.zBinSize);
        if(z<0) ++zbin;
        
        // std::cout<<"this triplet has z = "<<z<<" and will go to the bin "<<zbin<<std::endl;

        // extend the histogram with zeros, if needed 
        while(zbin > hist.size()-1)
        {
            hist.resize(2*hist.size(),0);
            // ACTS_DEBUG("new hist size = "<<hist.size());
        }

        ++hist.at(zbin);
    }

    return hist;
}


template <typename spacepoint_t>
float Acts::ZScanSeedVertexFinder<spacepoint_t>::findZPeak(const std::vector<int>& hist, int doSecondPeak) const
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
        float z=m_cfg.zBinSize*((h/2)+0.5);
        // increase/decrease the position by a half bin size, so it corresponds to the middle of the bin
        if(h%2)
        {
            // negative z value for odd bins
            z*=-1;
        }
        // make weighted average
        zpos+=z*hist.at(h);
        zsum+=hist.at(h);
    }

    if(doSecondPeak==0 || doSecondPeak==1)
    {
        if(!doSecondPeak) return zpos/zsum;

        std::cout<<"found one peak at "<<zpos/zsum<<std::endl;

        float minSecondPeakMagnitude=0.8;
        int maxh2=0;
        for(unsigned int h2=0;h2<hist.size();h2++)
        {
            if(maxh!=1)
            {
                // works for both maxh==0 and maxh>=2
                if(std::abs(maxh-(int)h2)<3) continue;
            }
            else
            {
                if(h2==3 || h2==0) continue;
            }

            if(hist.at(h2) > hist.at(maxh2)) maxh2=h2;
        }

        // look around the maximum for a better estimation of the peak
        std::vector<int> h2peak{maxh2-2,maxh2,maxh2+2};
        // fix for bins around zero
        if(maxh2==0) h2peak.at(0)=1;
        if(maxh2==1) h2peak.at(0)=0;

        float zsum2=0.,zpos2=0.;
        for(auto h2 : h2peak)
        {
            // position of the edge of the bin closer to zero
            float z2=m_cfg.zBinSize*((h2/2)+0.5);
            // increase/decrease the position by a half bin size, so it corresponds to the middle of the bin
            if(h2%2)
            {
                // negative z value for odd bins
                z2*=-1;
            }
            // make weighted average
            zpos2+=z2*hist.at(h2);
            zsum2+=hist.at(h2);
        }

        if(zsum2 > minSecondPeakMagnitude*zsum)
        {
            std::cout<<"found second peak at "<<zpos2/zsum2<<std::endl;
            return (zpos+zpos2)/(zsum+zsum2);
        }
        else
        {
            return zpos/zsum;
        }
    }

    if(doSecondPeak==2)
    {
        // find mean
        // hist.at(maxh)

        float zposM=0., zsumM=0.;
        for(unsigned int h=0;h<hist.size();h++)
        {
            if(hist.at(h)<0.1*hist.at(maxh)) continue;
            
            // increase the position by a half bin size, so it corresponds to the middle of the bin
            float z=m_cfg.zBinSize*((h/2)+0.5);
            if(h%2)
            {
                // negative z value for odd bins
                z*=-1;
            }   
            
            // make weighted average
            zposM+=z*hist.at(h);
            zsumM+=hist.at(h);
        }

        return zposM/zsumM;
    }

    return zpos/zsum;
}
