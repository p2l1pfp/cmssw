#ifndef __L1Trigger_L1THGCal_HGCalHistoClustering_h__
#define __L1Trigger_L1THGCal_HGCalHistoClustering_h__

#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalTriggerCell_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistogramCell_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/CentroidHelper.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringConfig_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalCluster_SA.h"

namespace l1thgcfirmware {

  class HGCalHistoClustering {
  public:
    HGCalHistoClustering(const l1thgcfirmware::ClusterAlgoConfig& config);
    ~HGCalHistoClustering() {}

    void runClustering(const l1thgcfirmware::HGCalTriggerCellSAPtrCollection& triggerCellsIn,
                       const l1thgcfirmware::HGCalHistogramCellSAPtrCollection& histogramIn,
                       l1thgcfirmware::HGCalTriggerCellSAShrPtrCollection& clusteredTriggerCellsOut,
                       l1thgcfirmware::CentroidHelperPtrCollection& readoutFlagsOut,
                       l1thgcfirmware::HGCalClusterSAPtrCollection& protoClustersOut) const;

  private:
    // Clustering
    void clusterizer(const l1thgcfirmware::HGCalTriggerCellSAPtrCollection& triggerCellsIn,
                     const l1thgcfirmware::HGCalHistogramCellSAPtrCollection& histogram,
                     l1thgcfirmware::HGCalTriggerCellSAShrPtrCollection& clusteredTriggerCells,
                     l1thgcfirmware::HGCalTriggerCellSAShrPtrCollection& unclusteredTriggerCells,
                     l1thgcfirmware::CentroidHelperPtrCollection& readoutFlags) const;
    void triggerCellToCluster(const l1thgcfirmware::HGCalTriggerCellSAShrPtrCollection& clusteredTriggerCells,
                              l1thgcfirmware::HGCalClusterSAPtrCollection& clustersOut) const;

    const l1thgcfirmware::ClusterAlgoConfig& config_;
  };
}  // namespace l1thgcfirmware

#endif