#ifndef __L1Trigger_L1THGCal_HGCalTCDistribution_h__
#define __L1Trigger_L1THGCal_HGCalTCDistribution_h__

#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalTriggerCell_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/DistServer.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringConfig_SA.h"

#include <vector>

namespace l1thgcfirmware {

  class HGCalTCDistribution {
  public:
    HGCalTCDistribution(const l1thgcfirmware::ClusterAlgoConfig& config);
    ~HGCalTCDistribution() {}

    void runTriggerCellDistribution(const l1thgcfirmware::HGCalTriggerCellSAPtrCollections& triggerCellsIn,
                                    l1thgcfirmware::HGCalTriggerCellSAPtrCollection& triggerCellsOut) const;

  private:
    // TC input step
    void triggerCellInput(const l1thgcfirmware::HGCalTriggerCellSAPtrCollections& inputTCs,
                          l1thgcfirmware::HGCalTriggerCellSAShrPtrCollection& outputTCs) const;

    // TC distribution steps
    void triggerCellDistribution0(const l1thgcfirmware::HGCalTriggerCellSAShrPtrCollection& triggerCellsIn) const;
    void triggerCellDistribution1(
        const l1thgcfirmware::HGCalTriggerCellSAShrPtrCollection& triggerCellsIn,
        l1thgcfirmware::HGCalTriggerCellSAShrPtrCollections& outTriggerCellDistributionGrid) const;
    void triggerCellDistribution2(
        const l1thgcfirmware::HGCalTriggerCellSAShrPtrCollections& inTriggerCellDistributionGrid,
        l1thgcfirmware::HGCalTriggerCellSAShrPtrCollection& triggerCellsOut,
        l1thgcfirmware::HGCalTriggerCellSAShrPtrCollections& outTriggerCellDistributionGrid) const;
    void triggerCellDistribution3(
        const l1thgcfirmware::HGCalTriggerCellSAShrPtrCollections& inTriggerCellDistributionGrid,
        l1thgcfirmware::HGCalTriggerCellSAShrPtrCollection& triggerCellsOut,
        l1thgcfirmware::HGCalTriggerCellSAShrPtrCollections& outTriggerCellDistributionGrid) const;
    void triggerCellDistribution4(const l1thgcfirmware::HGCalTriggerCellSAShrPtrCollection& triggerCellsIn) const;
    void triggerCellDistribution5(
        const l1thgcfirmware::HGCalTriggerCellSAShrPtrCollections& inTriggerCellDistributionGrid,
        l1thgcfirmware::HGCalTriggerCellSAPtrCollection& triggerCellsOut) const;

    // Useful functions
    void initializeTriggerCellDistGrid(l1thgcfirmware::HGCalTriggerCellSAShrPtrCollections& grid,
                                       unsigned int nX,
                                       unsigned int nY) const;

    void runDistServers(const l1thgcfirmware::HGCalTriggerCellSAShrPtrCollections& gridIn,
                        l1thgcfirmware::HGCalTriggerCellSAShrPtrCollections& gridOut,
                        l1thgcfirmware::HGCalTriggerCellSAShrPtrCollection& tcsOut,
                        unsigned int latency,
                        unsigned int nDistServers,
                        unsigned int nInputs,
                        unsigned int nOutputs,
                        unsigned int nInterleave,
                        bool setOutputGrid) const;

    const l1thgcfirmware::ClusterAlgoConfig& config_;
  };
}  // namespace l1thgcfirmware

#endif