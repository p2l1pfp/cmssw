#ifndef __L1Trigger_L1THGCal_HGCalHistoCluteringConfig_SA_h__
#define __L1Trigger_L1THGCal_HGCalHistoCluteringConfig_SA_h__

#include <map>
#include <vector>

namespace l1thgcfirmware {

  enum Step {
    Uninitialized = -1,
    Input,
    Dist0,
    Dist1,
    Dist2,
    Dist3,
    Dist4,
    Dist5,
    TcToHc,
    Hist,
    Smearing1D,
    NormArea,
    Smearing2D,
    Maxima1D,  // Not actually used currently
    Maxima2D,
    CalcAverage,
    Clusterizer,
    TriggerCellToCluster,
    ClusterSum
  };

  class ClusterAlgoConfig {
  public:
    ClusterAlgoConfig();

    void setParameters() {}

    void setSector(const unsigned int sector) { sector_ = sector; }
    unsigned int sector() const { return sector_; }

    void setZSide(const int zSide) { zSide_ = zSide; }
    int zSide() const { return zSide_; }

    void setStepLatencies(const std::vector<unsigned int>& latencies);
    unsigned int getStepLatency(const Step step) const { return stepLatency_.at(step); }
    unsigned int getLatencyUpToAndIncluding(const Step step) const;

    void setClusterizerOffset(const unsigned clusterizerOffset) { clusterizerOffset_ = clusterizerOffset; }
    unsigned int clusterizerOffset() const { return clusterizerOffset_; }

    void setCClocks(const unsigned cClocks) { cClocks_ = cClocks; }
    unsigned int cClocks() const { return cClocks_; }

    void setCInputs(const unsigned cInputs) { cInputs_ = cInputs; }
    unsigned int cInputs() const { return cInputs_; }

    void setCInputs2(const unsigned cInputs2) { cInputs2_ = cInputs2; }
    unsigned int cInputs2() const { return cInputs2_; }

    void setCInt(const unsigned cInt) { cInt_ = cInt; }
    unsigned int cInt() const { return cInt_; }

    void setCColumns(const unsigned cColumns) { cColumns_ = cColumns; }
    unsigned int cColumns() const { return cColumns_; }

    void setCRows(const unsigned cRows) { cRows_ = cRows; }
    unsigned int cRows() const { return cRows_; }

    void setROverZHistOffset(const unsigned rOverZHistOffset) { rOverZHistOffset_ = rOverZHistOffset; }
    unsigned int rOverZHistOffset() const { return rOverZHistOffset_; }

    void setROverZBinSize(const unsigned rOverZBinSize) { rOverZBinSize_ = rOverZBinSize; }
    unsigned int rOverZBinSize() const { return rOverZBinSize_; }

    void setMaxBinsSmearing1D(const unsigned maxBinsSmearing1D) { maxBinsSmearing1D_ = maxBinsSmearing1D; }
    unsigned maxBinsSmearing1D() const { return maxBinsSmearing1D_; }
    void setNBitsAreaNormLUT(const unsigned nBitsAreaNormLUT) { nBitsAreaNormLUT_ = nBitsAreaNormLUT; }
    unsigned nBitsAreaNormLUT() const { return nBitsAreaNormLUT_; }
    unsigned int kernelWidth(unsigned int iBin) const { return kernelWidths_.at(iBin); }
    unsigned int areaNormalization(unsigned int iBin) const { return areaNormalizations_.at(iBin); }

    void setROverZRange(const float rOverZRange) { rOverZRange_ = rOverZRange; }
    float rOverZRange() const { return rOverZRange_; }
    void setROverZNValues(const unsigned rOverZNValues) { rOverZNValues_ = rOverZNValues; }
    unsigned rOverZNValues() const { return rOverZNValues_; }
    void setPhiRange(const float phiRange) { phiRange_ = phiRange; }
    float phiRange() const { return phiRange_; }
    void setPhiNValues(const unsigned phiNValues) { phiNValues_ = phiNValues; }
    unsigned phiNValues() const { return phiNValues_; }
    void setPtDigiFactor(const float ptDigiFactor) { ptDigiFactor_ = ptDigiFactor; }
    float ptDigiFactor() const { return ptDigiFactor_; }

    void setMinClusterPtOut(const float pt) { minClusterPtOut_ = pt; }
    float minClusterPtOut() const { return minClusterPtOut_; }

    void setMaxClustersPerLink(const unsigned maxClustersPerLink) { maxClustersPerLink_ = maxClustersPerLink; }
    unsigned maxClustersPerLink() const { return maxClustersPerLink_; }
    void setNInputLinks(const unsigned nInputLinks) { nInputLinks_ = nInputLinks; }
    unsigned nInputLinks() const { return nInputLinks_; }

    void setN60Sectors(const unsigned n60Sectors) { n60Sectors_ = n60Sectors; }
    unsigned n60Sectors() const { return n60Sectors_; }
    void setNCoarsePhiDist1(const unsigned nCoarsePhiDist1) { nCoarsePhiDist1_ = nCoarsePhiDist1; }
    unsigned nCoarsePhiDist1() const { return nCoarsePhiDist1_; }
    void setNDistServers1(const unsigned nDistServers1) { nDistServers1_ = nDistServers1; }
    unsigned nDistServers1() const { return nDistServers1_; }
    void setDistServer1_nIn(const unsigned distServer1_nIn) { distServer1_nIn_ = distServer1_nIn; }
    unsigned distServer1_nIn() const { return distServer1_nIn_; }
    void setDistServer1_nOut(const unsigned distServer1_nOut) { distServer1_nOut_ = distServer1_nOut; }
    unsigned distServer1_nOut() const { return distServer1_nOut_; }
    void setDistServer1_nInterleave(const unsigned distServer1_nInterleave) {
      distServer1_nInterleave_ = distServer1_nInterleave;
    }
    unsigned distServer1_nInterleave() const { return distServer1_nInterleave_; }
    void setNCoarsePhiDist2(const unsigned nCoarsePhiDist2) { nCoarsePhiDist2_ = nCoarsePhiDist2; }
    unsigned nCoarsePhiDist2() const { return nCoarsePhiDist2_; }
    void setNDistServers2(const unsigned nDistServers2) { nDistServers2_ = nDistServers2; }
    unsigned nDistServers2() const { return nDistServers2_; }
    void setDistServer2_nIn(const unsigned distServer2_nIn) { distServer2_nIn_ = distServer2_nIn; }
    unsigned distServer2_nIn() const { return distServer2_nIn_; }
    void setDistServer2_nOut(const unsigned distServer2_nOut) { distServer2_nOut_ = distServer2_nOut; }
    unsigned distServer2_nOut() const { return distServer2_nOut_; }
    void setDistServer2_nInterleave(const unsigned distServer2_nInterleave) {
      distServer2_nInterleave_ = distServer2_nInterleave;
    }
    unsigned distServer2_nInterleave() const { return distServer2_nInterleave_; }

    void setClusterizerMagicTime(const unsigned clusterizerMagicTime) { clusterizerMagicTime_ = clusterizerMagicTime; }
    unsigned int clusterizerMagicTime() const { return clusterizerMagicTime_; }
    void setNFifos(const unsigned nFifos) { nFifos_ = nFifos; }
    unsigned int nFifos() const { return nFifos_; }
    void setNColumnsPerFifo(const unsigned nColumnsPerFifo) { nColumnsPerFifo_ = nColumnsPerFifo; }
    unsigned int nColumnsPerFifo() const { return nColumnsPerFifo_; }
    void setFirstSeedBin(const unsigned firstSeedBin) { firstSeedBin_ = firstSeedBin; }
    unsigned int firstSeedBin() const { return firstSeedBin_; }
    void setNColumnFifoVeto(const unsigned nColumnFifoVeto) { nColumnFifoVeto_ = nColumnFifoVeto; }
    unsigned int nColumnFifoVeto() const { return nColumnFifoVeto_; }
    void setDeltaR2Cut(const unsigned deltaR2Cut) { deltaR2Cut_ = deltaR2Cut; }
    unsigned int deltaR2Cut() const { return deltaR2Cut_; }
    void setNColumnsForClustering(const unsigned nColumnsForClustering) {
      nColumnsForClustering_ = nColumnsForClustering;
    }
    unsigned int nColumnsForClustering() const { return nColumnsForClustering_; }
    void setNRowsForClustering(const unsigned nRowsForClustering) { nRowsForClustering_ = nRowsForClustering; }
    unsigned int nRowsForClustering() const { return nRowsForClustering_; }

    void setThresholdParams(const unsigned int a, const unsigned int b, const int c) {
      thresholdMaximaParam_a_ = a;
      thresholdMaximaParam_b_ = b;
      thresholdMaximaParam_c_ = c;
    }
    unsigned int thresholdMaxima(unsigned int iBin) const { return thresholdMaximaConstants_.at(iBin); }

    void setNBinsCosLUT(const unsigned int nBinsCosLUT) { nBinsCosLUT_ = nBinsCosLUT; }
    unsigned int nBinsCosLUT() const { return nBinsCosLUT_; }
    void setNBitsCosLUT(const unsigned int nBitsCosLUT) { nBitsCosLUT_ = nBitsCosLUT; }
    unsigned int nBitsCosLUT() const { return nBitsCosLUT_; }

    unsigned int cosLUT(unsigned int iBin) const { return cosLUT_.at(iBin); }

    void setDepths(const std::vector<unsigned int>& depths) {
      depths_.clear();
      depths_.reserve(depths.size());
      for (const auto& depth : depths)
        depths_.push_back(depth);
    }
    const std::vector<unsigned int>& depths() const { return depths_; }
    unsigned int depth(unsigned int iLayer) const { return depths_.at(iLayer); }

    void setTriggerLayers(const std::vector<unsigned int>& triggerLayers) { triggerLayers_ = triggerLayers; }
    const std::vector<unsigned int>& triggerLayers() const { return triggerLayers_; }
    unsigned int triggerLayer(unsigned int iLayer) const { return triggerLayers_.at(iLayer); }

    void setLayerWeights_E(const std::vector<unsigned int> layerWeights_E) {
      layerWeights_E_.clear();
      layerWeights_E_.reserve(layerWeights_E.size());
      for (const auto& weight : layerWeights_E)
        layerWeights_E_.push_back(weight);
    }
    const std::vector<unsigned int>& layerWeights_E() const { return layerWeights_E_; }
    unsigned int layerWeight_E(unsigned int iTriggerLayer) const { return layerWeights_E_.at(iTriggerLayer); }

    void setLayerWeights_E_EM(const std::vector<unsigned int> layerWeights_E_EM) {
      layerWeights_E_EM_.clear();
      layerWeights_E_EM_.reserve(layerWeights_E_EM.size());
      for (const auto& weight : layerWeights_E_EM)
        layerWeights_E_EM_.push_back(weight);
    }
    const std::vector<unsigned int>& layerWeights_E_EM() const { return layerWeights_E_EM_; }
    unsigned int layerWeight_E_EM(unsigned int iTriggerLayer) const { return layerWeights_E_EM_.at(iTriggerLayer); }

    void setLayerWeights_E_EM_core(const std::vector<unsigned int> layerWeights_E_EM_core) {
      layerWeights_E_EM_core_.clear();
      layerWeights_E_EM_core_.reserve(layerWeights_E_EM_core.size());
      for (const auto& weight : layerWeights_E_EM_core)
        layerWeights_E_EM_core_.push_back(weight);
    }
    const std::vector<unsigned int>& layerWeights_E_EM_core() const { return layerWeights_E_EM_core_; }
    unsigned int layerWeight_E_EM_core(unsigned int iTriggerLayer) const {
      return layerWeights_E_EM_core_.at(iTriggerLayer);
    }

    void setLayerWeights_E_H_early(const std::vector<unsigned int> layerWeights_E_H_early) {
      layerWeights_E_H_early_.clear();
      layerWeights_E_H_early_.reserve(layerWeights_E_H_early.size());
      for (const auto& weight : layerWeights_E_H_early)
        layerWeights_E_H_early_.push_back(weight);
    }
    const std::vector<unsigned int>& layerWeights_E_H_early() const { return layerWeights_E_H_early_; }
    unsigned int layerWeight_E_H_early(unsigned int iTriggerLayer) const {
      return layerWeights_E_H_early_.at(iTriggerLayer);
    }

    void setCorrection(const unsigned correction) { correction_ = correction; }
    unsigned int correction() const { return correction_; }

    void setSaturation(const unsigned saturation) { saturation_ = saturation; }
    unsigned int saturation() const { return saturation_; }

    void setNTriggerLayers(const unsigned n) { nTriggerLayers_ = n; }
    unsigned int nTriggerLayers() const { return nTriggerLayers_; }

    void initializeLUTs();

    void printConfiguration();  // For debugging

  private:
    void initializeSmearingKernelConstants(unsigned int bins, unsigned int offset, unsigned int height);
    void initializeThresholdMaximaConstants(unsigned int bins, unsigned int a, unsigned int b, int c);
    void initializeCosLUT();

    unsigned int histogramOffset_;
    unsigned int clusterizerOffset_;
    unsigned int cClocks_;
    unsigned int cInputs_;
    unsigned int cInputs2_;  // Better name for variable?
    unsigned int cInt_;
    unsigned int cColumns_;
    unsigned int cRows_;
    unsigned int rOverZHistOffset_;
    unsigned int rOverZBinSize_;

    // Digitisation parameters
    float rOverZRange_;
    unsigned rOverZNValues_;
    float phiRange_;
    unsigned phiNValues_;
    float ptDigiFactor_;

    // Selection on output clusters
    float minClusterPtOut_;

    // Input link params
    unsigned int maxClustersPerLink_;
    unsigned int nInputLinks_;

    // TC distribution parameters
    unsigned int n60Sectors_;
    unsigned int nCoarsePhiDist1_;
    unsigned int nDistServers1_;
    unsigned int distServer1_nIn_;
    unsigned int distServer1_nOut_;
    unsigned int distServer1_nInterleave_;
    unsigned int nCoarsePhiDist2_;
    unsigned int nDistServers2_;
    unsigned int distServer2_nIn_;
    unsigned int distServer2_nOut_;
    unsigned int distServer2_nInterleave_;

    // Smearing and normalization parameters
    unsigned int maxBinsSmearing1D_;
    unsigned int nBitsAreaNormLUT_;
    std::vector<unsigned int> kernelWidths_;
    std::vector<unsigned int> areaNormalizations_;

    // Threshold maxima parameters
    unsigned int thresholdMaximaParam_a_;
    unsigned int thresholdMaximaParam_b_;
    int thresholdMaximaParam_c_;
    std::vector<int> thresholdMaximaConstants_;

    // Clusterizer parameters
    unsigned int nBinsCosLUT_;
    unsigned int nBitsCosLUT_;
    std::vector<unsigned int> cosLUT_;
    unsigned int nFifos_;
    unsigned int nColumnsPerFifo_;
    unsigned int firstSeedBin_;
    unsigned int nColumnFifoVeto_;
    unsigned int deltaR2Cut_;
    unsigned int nColumnsForClustering_;
    unsigned int nRowsForClustering_;
    unsigned int clusterizerMagicTime_;

    std::vector<unsigned int> stepLatency_;

    // Parameters for triggerCellToCluster
    std::vector<unsigned int> depths_;
    std::vector<unsigned int> triggerLayers_;
    std::vector<unsigned int> layerWeights_E_;
    std::vector<unsigned int> layerWeights_E_EM_;
    std::vector<unsigned int> layerWeights_E_EM_core_;
    std::vector<unsigned int> layerWeights_E_H_early_;
    unsigned int correction_;
    unsigned int saturation_;

    // Trigger geometry info
    unsigned int nTriggerLayers_;

    unsigned int sector_;
    int zSide_;
  };

}  // namespace l1thgcfirmware

#endif
