#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusterProperties_SA.h"

#include <cmath>
#include <algorithm>

using namespace std;
using namespace l1thgcfirmware;

HGCalHistoClusterProperties::HGCalHistoClusterProperties(const ClusterAlgoConfig& config) : config_(config) {}

void HGCalHistoClusterProperties::runClusterProperties(
    const l1thgcfirmware::HGCalClusterSAPtrCollection& protoClustersIn,
    const CentroidHelperPtrCollection& readoutFlags,
    HGCalClusterSAPtrCollection& clustersOut) const {
  // Cluster properties
  HGCalClusterSAPtrCollection clusterAccumulation;
  clusterSum(protoClustersIn, readoutFlags, clusterAccumulation, clustersOut);
  clusterProperties(clustersOut);
}

// Accumulates/combines inputs cluster objects (each corresponding to one TC belonging to a cluster) into clusters  (one per cluster made up of TCs)
void HGCalHistoClusterProperties::clusterSum(const HGCalClusterSAPtrCollection& protoClusters,
                                             const CentroidHelperPtrCollection& readoutFlags,
                                             HGCalClusterSAPtrCollection& clusterAccumulation,
                                             HGCalClusterSAPtrCollection& clusterSums) const {
  HGCalClusterSAPtrCollections protoClustersPerColumn(config_.cColumns());
  vector<unsigned int> clock(config_.cColumns(), 0);
  for (const auto& protoCluster : protoClusters) {
    auto index = protoCluster->index();
    // Do we have to make a copy of protoCluster here?
    protoClustersPerColumn.at(index).push_back(make_unique<HGCalCluster>(*protoCluster));
  }

  map<unsigned int, HGCalClusterSAPtr> sums;

  for (const auto& flag : readoutFlags) {
    auto accumulator = make_unique<HGCalCluster>(0, 0, true, true);
    const unsigned stepLatency = 23;
    flag->setClock(flag->clock() + stepLatency);

    for (const auto& protoCluster : protoClustersPerColumn.at(flag->index())) {
      if (protoCluster->clock() <= clock.at(flag->index()))
        continue;
      if (protoCluster->clock() > flag->clock())
        continue;
      *accumulator += *protoCluster;
    }

    clock.at(flag->index()) = flag->clock();
    accumulator->setClock(flag->clock());
    accumulator->setIndex(flag->index());
    accumulator->setDataValid(true);

    if (sums.find(flag->clock()) == sums.end()) {
      const unsigned stepLatency = 7;
      auto sum = make_unique<HGCalCluster>(flag->clock() + stepLatency, 0, true, true);
      sums[flag->clock()] = move(sum);
    }

    *(sums.at(flag->clock())) += *accumulator;

    clusterAccumulation.push_back(move(accumulator));
  }

  for (auto& sum : sums) {
    clusterSums.push_back(move(sum.second));
  }
}

// Calculates properties of clusters from accumulated quantities
void HGCalHistoClusterProperties::clusterProperties(HGCalClusterSAPtrCollection& clusterSums) const {
  unsigned int nTCs = 0;
  for (auto& c : clusterSums) {
    if (c->n_tc_w() == 0)
      continue;
    std::pair<unsigned int, unsigned int> sigmaEnergy = sigma_energy(c->n_tc_w(), c->w2(), c->w());
    c->set_sigma_e_quotient(sigmaEnergy.first);
    c->set_sigma_e_fraction(sigmaEnergy.second);
    std::pair<unsigned int, unsigned int> Mean_z = mean_coordinate(c->wz(), c->w());
    c->set_mean_z_quotient(Mean_z.first);
    c->set_mean_z_fraction(Mean_z.second);
    std::pair<unsigned int, unsigned int> Mean_phi = mean_coordinate(c->wphi(), c->w());
    c->set_mean_phi_quotient(Mean_phi.first);
    c->set_mean_phi_fraction(Mean_phi.second);
    std::pair<unsigned int, unsigned int> Mean_eta = mean_coordinate(c->weta(), c->w());
    c->set_mean_eta_quotient(Mean_eta.first);
    c->set_mean_eta_fraction(Mean_eta.second);
    std::pair<unsigned int, unsigned int> Mean_roz = mean_coordinate(c->wroz(), c->w());
    c->set_mean_roz_quotient(Mean_roz.first);
    c->set_mean_roz_fraction(Mean_roz.second);
    std::pair<unsigned int, unsigned int> Sigma_z = sigma_coordinate(c->w(), c->wz2(), c->wz());
    c->set_sigma_z_quotient(Sigma_z.first);
    c->set_sigma_z_fraction(Sigma_z.second);
    std::pair<unsigned int, unsigned int> Sigma_phi = sigma_coordinate(c->w(), c->wphi2(), c->wphi());
    c->set_sigma_phi_quotient(Sigma_phi.first);
    c->set_sigma_phi_fraction(Sigma_phi.second);
    std::pair<unsigned int, unsigned int> Sigma_eta = sigma_coordinate(c->w(), c->weta2(), c->weta());
    c->set_sigma_eta_quotient(Sigma_eta.first);
    c->set_sigma_eta_fraction(Sigma_eta.second);
    std::pair<unsigned int, unsigned int> Sigma_roz = sigma_coordinate(c->w(), c->wroz2(), c->wroz());
    c->set_sigma_roz_quotient(Sigma_roz.first);
    c->set_sigma_roz_fraction(Sigma_roz.second);
    std::vector<int> layeroutput = showerLengthProperties(c->layerbits());
    c->set_firstLayer(layeroutput[0]);
    c->set_lastLayer(layeroutput[1]);
    c->set_showerLen(layeroutput[2]);
    c->set_coreShowerLen(layeroutput[3]);
    std::pair<unsigned int, unsigned int> e_em_over_e = energy_ratio(c->e_em(), c->e());
    c->set_e_em_over_e_quotient(e_em_over_e.first);
    c->set_e_em_over_e_fraction(e_em_over_e.second);
    std::pair<unsigned int, unsigned int> e_em_core_over_e_em = energy_ratio(c->e_em_core(), c->e());
    c->set_e_em_core_over_e_em_quotient(e_em_core_over_e_em.first);
    c->set_e_em_core_over_e_em_fraction(e_em_core_over_e_em.second);
    std::pair<unsigned int, unsigned int> e_h_early_over_e = energy_ratio(c->e_h_early(), c->e());
    c->set_e_h_early_over_e_quotient(e_h_early_over_e.first);
    c->set_e_h_early_over_e_fraction(e_h_early_over_e.second);

    nTCs += c->n_tc();
  }
}

std::pair<unsigned int, unsigned int> HGCalHistoClusterProperties::sigma_energy(unsigned int N_TC_W,
                                                                                unsigned long int Sum_W2,
                                                                                unsigned int Sum_W) const {
  unsigned long int N = N_TC_W * Sum_W2 - pow(Sum_W, 2);
  unsigned long int D = pow(N_TC_W, 2);
  if (D == 0) {
    return {0, 0};
  }
  double intpart;
  const unsigned shift = 2;  // Shift by one bit, pow(2,1)
  double frac = modf(sqrt(N / D), &intpart) * shift;
  return {(unsigned int)intpart, (unsigned int)frac};
}

std::pair<unsigned int, unsigned int> HGCalHistoClusterProperties::mean_coordinate(unsigned int Sum_Wc,
                                                                                   unsigned int Sum_W) const {
  if (Sum_W == 0) {
    return {0, 0};
  }
  double intpart;
  const unsigned shift = 4;  // Shift by one bit, pow(2,2)
  double frac = modf((double)Sum_Wc / Sum_W, &intpart) * shift;
  return {(unsigned int)intpart, (unsigned int)frac};
}

std::pair<unsigned int, unsigned int> HGCalHistoClusterProperties::sigma_coordinate(unsigned int Sum_W,
                                                                                    unsigned long int Sum_Wc2,
                                                                                    unsigned int Sum_Wc) const {
  unsigned long int N = Sum_W * Sum_Wc2 - pow(Sum_Wc, 2);
  unsigned long int D = pow(Sum_W, 2);
  if (D == 0) {
    return {0, 0};
  }
  double intpart;
  const unsigned shift = 2;  // Shift by one bit, pow(2,1)
  double frac = modf((double)sqrt(N / D), &intpart) * shift;
  return {(unsigned int)intpart, (unsigned int)frac};
}

std::pair<unsigned int, unsigned int> HGCalHistoClusterProperties::energy_ratio(unsigned int e_N,
                                                                                unsigned int e_D) const {
  if (e_D == 0) {
    return {0, 0};
  } else {
    double intpart;
    const unsigned shift = 256;  // Shift by eight bit, pow(2,8)
    double frac = modf((double)e_N / e_D, &intpart) * shift;
    return {(unsigned int)intpart, (unsigned int)frac};
  }
}

std::vector<int> HGCalHistoClusterProperties::showerLengthProperties(unsigned long int layerBits) const {
  int counter = 0;
  int firstLayer = 0;
  bool firstLayerFound = false;
  int lastLayer = 0;
  std::vector<int> layerBits_array;

  for (unsigned int idx = 0; idx < config_.nTriggerLayers(); idx++) {
    if ((layerBits & (1L << (config_.nTriggerLayers() - 1 - idx))) >= 1L) {
      if (!firstLayerFound) {
        firstLayer = idx + 1;
        firstLayerFound = true;
      }
      lastLayer = idx + 1;
      counter += 1;
    } else {
      layerBits_array.push_back(counter);
      counter = 0;
    }
  }
  int showerLen = lastLayer - firstLayer + 1;
  int coreShowerLen = config_.nTriggerLayers();
  if (!layerBits_array.empty()) {
    coreShowerLen = *std::max_element(layerBits_array.begin(), layerBits_array.end());
  }

  return {firstLayer, lastLayer, showerLen, coreShowerLen};
}
