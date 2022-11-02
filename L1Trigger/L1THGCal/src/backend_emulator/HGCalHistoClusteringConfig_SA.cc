#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringConfig_SA.h"
#include <cmath>
#include <iostream>
using namespace std;
using namespace l1thgcfirmware;
ClusterAlgoConfig::ClusterAlgoConfig()
    : clusterizerOffset_(0),
      cClocks_(0),
      cInputs_(0),
      cInputs2_(0),
      cInt_(0),
      cColumns_(0),
      cRows_(0),
      rOverZHistOffset_(0),
      rOverZBinSize_(0),
      kernelWidths_(),
      areaNormalizations_(),
      thresholdMaximaConstants_(),
      nBinsCosLUT_(0),
      cosLUT_(),
      clusterizerMagicTime_(0),
      stepLatency_(),
      depths_(),
      triggerLayers_(),
      layerWeights_E_(),
      layerWeights_E_EM_(),
      layerWeights_E_EM_core_(),
      layerWeights_E_H_early_(),
      correction_(),
      saturation_() {}

void ClusterAlgoConfig::setStepLatencies(const std::vector<unsigned int>& latencies) {
  stepLatency_.clear();
  stepLatency_.reserve(latencies.size());
  for (unsigned int iStep = 0; iStep < latencies.size(); ++iStep) {
    stepLatency_.push_back(latencies.at(iStep));
  }
}

unsigned int ClusterAlgoConfig::getLatencyUpToAndIncluding(const Step step) const {
  unsigned int latency = 0;
  for (int iStep = 0; iStep <= step; ++iStep)
    latency += getStepLatency(Step(iStep));
  return latency;
}

void ClusterAlgoConfig::initializeLUTs() {
  initializeSmearingKernelConstants(cRows_, rOverZHistOffset_, rOverZBinSize_);
  initializeThresholdMaximaConstants(cRows_, thresholdMaximaParam_a_, thresholdMaximaParam_b_, thresholdMaximaParam_c_);
  initializeCosLUT();
}

void ClusterAlgoConfig::initializeSmearingKernelConstants(unsigned int bins, unsigned int offset, unsigned int height) {
  const unsigned int lWidth0 = offset + (0.5 * height);
  const unsigned int lTarget = int((maxBinsSmearing1D_ + 0.5) * lWidth0 - 0.5);
  for (unsigned int iBin = 0; iBin < bins; ++iBin) {
    unsigned int lCentre = lWidth0 + (height * iBin);
    const unsigned int lBins = int(round(float(lTarget) / lCentre));

    kernelWidths_.push_back(lBins);

    lCentre *= lBins;

    const unsigned int lRatio = int(round(float(lTarget) / lCentre * (0x1 << nBitsAreaNormLUT_)));

    areaNormalizations_.push_back(lRatio);
  }
}

void ClusterAlgoConfig::initializeThresholdMaximaConstants(unsigned int bins, unsigned int a, unsigned int b, int c) {
  for (unsigned int iBin = 0; iBin < bins; ++iBin) {
    int threshold = a + b * iBin + c * iBin * iBin;
    thresholdMaximaConstants_.push_back(threshold);
  }
}

void ClusterAlgoConfig::initializeCosLUT() {
  for (unsigned int iBin = 0; iBin < nBinsCosLUT_ + 1; ++iBin) {
    unsigned int cosBin = round((0x1 << nBitsCosLUT_) * (1 - cos(iBin * M_PI / phiNValues_)));
    cosLUT_.push_back(cosBin);
  }
}
