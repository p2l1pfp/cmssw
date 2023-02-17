#include "L1Trigger/L1THGCalUtilities/interface/patternFiles/codecs_clusters.h"
#include <iostream>

namespace l1thgcfirmware {
    std::vector<HGCalCluster_HW> decodeClusters(const std::array<std::vector<ap_uint<64>>, 4>& frames) {

        std::vector<HGCalCluster_HW> clusters;
        unsigned nFrames = frames.at(0).size();
        for ( unsigned int iFrame=0; iFrame < nFrames; ++iFrame ) {
            if ( iFrame < 31 ) continue; // Skip header and towers

            ClusterWords clusterWords;
            for ( unsigned int i=0; i<nWordsPerCluster; ++i ) {
                clusterWords[i] = frames.at(i).at(iFrame);
            }
            HGCalCluster_HW cluster = HGCalCluster_HW::unpack( clusterWords );
            // std::cout << "Got a cluster on frame : " << iFrame << " " << cluster.e << std::endl;
            if (cluster.e > 0 ) clusters.emplace_back( cluster );
        }
        return clusters;
    }
}