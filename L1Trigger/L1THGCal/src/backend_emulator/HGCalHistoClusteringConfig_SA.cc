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
      saturation_(),
      rozToEtaLUT_() {}

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
  initializeRoZToEtaLUT();
  initializeSigmaRoZToSigmaEtaLUT();
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

void ClusterAlgoConfig::initializeThresholdMaximaConstants(unsigned int bins, unsigned int a, int b, int c) {
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

void ClusterAlgoConfig::initializeRoZToEtaLUT() {
  const float eta_min = 320 * M_PI/720; // ~1.4
  const float eta_max = 687 * M_PI/720; // ~3.0
  const float eta_LSB = M_PI/720;

  const float roz_min = 1/(sinh(eta_max)); // coming from L1T eta_max
  const float roz_max = 1/(sinh(eta_min)); // coming from L1T eta_min
  const float roz_LSB = (roz_max-roz_min)/pow(2,10); // coming from L1T eta_LSB; CP block needs to calculate mean(r/z) with at least 10b if we are to meet L1T eta LSB

  unsigned nEntries = 0x1<<10;
  rozToEtaLUT_.reserve(nEntries);
  for ( unsigned i = 0; i<nEntries; ++i ) {
    float roz = i * roz_LSB + roz_min;
    float eta = asinh(1./roz);
    unsigned eta_int_local = int( (eta - eta_min) / eta_LSB );
    unsigned eta_int_global = int( eta_int_local + eta_min / eta_LSB );
    rozToEtaLUT_.push_back(eta_int_global);
  }
}

void ClusterAlgoConfig::initializeSigmaRoZToSigmaEtaLUT() {

  // const float eta_min = 320 * M_PI/720; // ~1.4
  // const float eta_max = 687 * M_PI/720; // ~3.0
  // const float eta_LSB = M_PI/720;

  const float roz_min = 0.078953360909371;
  const float roz_max = 0.4869046060590794;
  const unsigned roz_nbits = 6;
  const float roz_LSB = (roz_max-roz_min)/pow(2,roz_nbits);

  // dEta/dMu LUT
  unsigned nEntries_dEta_dMu_LUT = 0x1<<roz_nbits;
  vector<double> dEta_dMu_LUT;
  dEta_dMu_LUT.reserve(nEntries_dEta_dMu_LUT);
  // std::cout << "dEta/dMu LUT" << std::endl;
  for ( unsigned i = 0; i<nEntries_dEta_dMu_LUT; ++i ) {
    float roz = i * roz_LSB + roz_min;
    float dEta_dMu = 1. / ( sqrt( pow(roz,-2)+1 ) * roz * roz );
    // std::cout << i << " " << roz_min << " " << roz_LSB << " " << roz << " " << dEta_dMu << std::endl;
    dEta_dMu_LUT.push_back(dEta_dMu);
  }

  // Sigma RoZ LUT
  const float sigmaRoz_min = 0;
  const float sigmaRoz_max = 0.0282874636266928;
  const unsigned sigmaRoz_nbits = 6;
  const float sigma_roz_LSB = (sigmaRoz_max-sigmaRoz_min)/pow(2,sigmaRoz_nbits);
  unsigned nEntries_sigmaRoz_LUT = 0x1<<sigmaRoz_nbits;
  vector<double> sigmaRoz_LUT;
  sigmaRoz_LUT.reserve(nEntries_sigmaRoz_LUT);
  for ( unsigned i = 0; i<nEntries_sigmaRoz_LUT; ++i ) {
    float sigma_roz = i * sigma_roz_LSB + sigmaRoz_min;
    sigmaRoz_LUT.push_back(sigma_roz);
  }

  // Final roz, sigma roz LUT
  const float sigmaEta_min = 0;
  const float sigmaEta_max = 0.1459855741275447;
  const unsigned sigmaEta_nbits = 5;
  const float sigma_eta_LSB = (sigmaEta_max-sigmaEta_min)/pow(2,sigmaEta_nbits);
  unsigned nEntries_roz_sigmaRoz_LUT = 0x1<<(roz_nbits+sigmaRoz_nbits);
  // vector<unsigned> roz_SigmaRoz_LUT(nEntries_roz_sigmaRoz_LUT);
  sigmaRozToSigmaEtaLUT_.reserve(nEntries_roz_sigmaRoz_LUT);
  for ( unsigned mu=0; mu<nEntries_dEta_dMu_LUT; ++mu) {
    double dEta_dMu = dEta_dMu_LUT.at(mu);
    for ( unsigned sigma=0; sigma<nEntries_sigmaRoz_LUT; ++sigma) {
      double sigmaRoz = sigmaRoz_LUT.at(sigma);

      float sigmaEta = dEta_dMu * sigmaRoz / sigma_eta_LSB;
      if ( sigmaEta > ( 0x1<<(sigmaEta_nbits-1) ) ) sigmaEta = 0x1<<(sigmaEta_nbits-1);
      sigmaRozToSigmaEtaLUT_.push_back(sigmaEta);
      // std::cout << mu << " " << sigma << " " << dEta_dMu << " " << sigmaRoz << " " << sigma_eta_LSB << " " << sigmaEta << std::endl;
    }
  }
}