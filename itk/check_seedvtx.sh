#!/bin/bash


# numParticles="100"
# particle="13"
# minPT="0.1 0.5 1.0"
# vtxX="0.0 0.5 1.0 1.5 2.0"
numParticles="426" #0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33"
# particle="1000822080"
part="1000822080"
# minPT="0.5"
pT="0.5"
vtxXlist="0.0"
vtxYlist="0.0"
# vtxX="0.5"
# vtxY="0.0"
# vtxZ="-180.89 -130.17 -80.31 -40.72 -10.66 20.48 50.51 90.03 140.26"
# vtxZlist="140.26"
vtxZlist="0.0"

cp -v ../Examples/Algorithms/Vertexing/src/SingleSeedVertexFinderAlgorithm.cpp SingleSeedVertexFinderAlgorithm.cpp_version31
cp -v ../Core/include/Acts/Vertexing/SingleSeedVertexFinder.ipp SingleSeedVertexFinder.ipp_version31


for numPart in `echo $numParticles`
    do
    for vtxX in `echo $vtxXlist`
        do
        for vtxY in `echo $vtxYlist`
            do
            for vtxZ in `echo $vtxZlist`
                do
                python3 full_chain_itk.py $numPart $part $pT $vtxX $vtxY $vtxZ
                mv -v itk_output_Hough/performance_seedvertexing.root performance_seedvertexing31_${numPart}_${part}_${pT}_${vtxX}_${vtxY}_${vtxZ}_Hough.root 
            done
        done
    done
done

# f reco track-associated truth particles in event : 9
# 16:11:21    RootVertexPe   INFO      Total number of reco track-associated truth primary vertices : 1
# 16:11:22    ZScanSeedVer   INFO      Found 1 vertices in the event in 1.48928 seconds.
# 16:11:22    ZScanSeedVer   INFO      Found vertex at z = 12.5306mm
# 16:11:22    ZScanSeedVer   INFO      Found 1 vertices in the event in 1.61118 seconds.
# 16:11:22    ZScanSeedVer   INFO      Found vertex at z = 14.7521mm
# 16:11:22    ZScanSeedVer   INFO      Found 1 vertices in the event in 2.00455 seconds.
# 16:11:22    ZScanSeedVer   INFO      Found vertex at z = 13.4082mm
# python3: /usr/include/boost/move/algo/adaptive_sort.hpp:452: bool boost::movelib::detail_adaptive::adaptive_sort_build_params(RandIt, Unsigned, Compare, Unsigned&, Unsigned&, Unsigned&, Unsigned&, XBuf&) [with RandIt = boost::container::dtl::pair<ActsFatras::Barcode, unsigned int>*; Compare = boost::container::dtl::flat_tree_value_compare<std::less<ActsFatras::Barcode>, boost::container::dtl::pair<ActsFatras::Barcode, unsigned int>, boost::container::dtl::select1st<ActsFatras::Barcode> >; Unsigned = long unsigned int; XBuf = boost::movelib::adaptive_xbuf<boost::container::dtl::pair<ActsFatras::Barcode, unsigned int>, boost::container::dtl::pair<ActsFatras::Barcode, unsigned int>*, long unsigned int>]: Assertion `collected < (n_min_ideal_keys+l_intbuf)' failed.
# ./check_seedvtx.sh: line 21: 15640 Aborted                 (core dumped) python3 full_chain_itk.py $numPart $part $pT $vtx $vtxY $vtxZ
# 16:11:23    buildITkGeom   INFO      Adding material from /ATLAS/pbalek/acts/acts/itk/acts-itk/itk-hgtd/material-maps-ITk-HGTD.json


