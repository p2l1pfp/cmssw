#ifndef __L1Trigger_L1THGCal_HGCalHistoSeeding_h__
#define __L1Trigger_L1THGCal_HGCalHistoSeeding_h__

#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalTriggerCell_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistogramCell_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringConfig_SA.h"

namespace l1thgcfirmware {

  class HGCalHistoSeeding {
  public:
    HGCalHistoSeeding(const l1thgcfirmware::ClusterAlgoConfig& config);
    ~HGCalHistoSeeding() {}

    void runSeeding(const l1thgcfirmware::HGCalTriggerCellSAPtrCollection& triggerCellsIn,
                    l1thgcfirmware::HGCalHistogramCellSAPtrCollection& histogramOut) const;

  private:
    // Histogram steps
    void triggerCellToHistogramCell(const l1thgcfirmware::HGCalTriggerCellSAPtrCollection& triggerCellsIn,
                                    l1thgcfirmware::HGCalHistogramCellSAPtrCollection& histogramOut) const;
    void makeHistogram(const l1thgcfirmware::HGCalHistogramCellSAPtrCollection& histogramCells,
                       l1thgcfirmware::HGCalHistogramCellSAPtrCollection& histogramOut) const;

    // Smearing steps
    void smearHistogram1D(l1thgcfirmware::HGCalHistogramCellSAPtrCollection& histogram) const;
    void normalizeArea(l1thgcfirmware::HGCalHistogramCellSAPtrCollection& histogram) const;
    void smearHistogram2D(l1thgcfirmware::HGCalHistogramCellSAPtrCollection& histogram) const;

    // Maxima finding
    void thresholdMaximaFinder(l1thgcfirmware::HGCalHistogramCellSAPtrCollection& histogram) const;
    void localMaximaFinder(l1thgcfirmware::HGCalHistogramCellSAPtrCollection& histogram) const;
    void calculateAveragePosition(l1thgcfirmware::HGCalHistogramCellSAPtrCollection& histogram) const;

    const l1thgcfirmware::ClusterAlgoConfig& config_;
  };
}  // namespace l1thgcfirmware

#endif