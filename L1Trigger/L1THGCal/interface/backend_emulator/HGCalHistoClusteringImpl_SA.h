#ifndef __L1Trigger_L1THGCal_HGCalHistoClusteringImplSA_h__
#define __L1Trigger_L1THGCal_HGCalHistoClusteringImplSA_h__

#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalTriggerCell_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistogramCell_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/CentroidHelper.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalCluster_SA.h"

#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalTCDistribution_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoSeeding_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClustering_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusterProperties_SA.h"

#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringConfig_SA.h"

namespace l1thgcfirmware {

  class HGCalHistoClusteringImplSA {
  public:
    HGCalHistoClusteringImplSA(const l1thgcfirmware::ClusterAlgoConfig& config);
    ~HGCalHistoClusteringImplSA() {}

    void runAlgorithm(const HGCalTriggerCellSAPtrCollections& inputs,
                      HGCalTriggerCellSAShrPtrCollection& clusteredTCs,
                      HGCalClusterSAPtrCollection& clusterSums) const;

  private:
    const l1thgcfirmware::ClusterAlgoConfig& config_;
    l1thgcfirmware::HGCalTCDistribution tcDistribution_;
    l1thgcfirmware::HGCalHistoSeeding seeding_;
    l1thgcfirmware::HGCalHistoClustering clustering_;
    l1thgcfirmware::HGCalHistoClusterProperties clusterProperties_;
  };

}  // namespace l1thgcfirmware

#endif
