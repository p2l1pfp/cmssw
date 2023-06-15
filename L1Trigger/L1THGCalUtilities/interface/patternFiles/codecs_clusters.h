#ifndef L1Trigger_L1TriggerUtilities_patternFiles_codecs_clusters_h
#define L1Trigger_L1TriggerUtilities_patternFiles_codecs_clusters_h

#include <array>
#include <vector>

#include "ap_int.h"

#include "DataFormats/Common/interface/View.h"
#include "L1Trigger/DemonstratorTools/interface/BoardData.h"
#include "DataFormats/L1THGCal/interface/HGCalCluster_HW.h"

namespace l1thgcfirmware {
    std::vector<HGCalCluster_HW> decodeClusters(const std::array<std::vector<ap_uint<64>>, 4> &);
}

#endif