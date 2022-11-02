#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringImpl_SA.h"

using namespace std;
using namespace l1thgcfirmware;
HGCalHistoClusteringImplSA::HGCalHistoClusteringImplSA(const ClusterAlgoConfig& config)
    : config_(config), tcDistribution_(config), seeding_(config), clustering_(config), clusterProperties_(config) {}

void HGCalHistoClusteringImplSA::runAlgorithm(const HGCalTriggerCellSAPtrCollections& inputs,
                                              HGCalTriggerCellSAShrPtrCollection& clusteredTCs,
                                              HGCalClusterSAPtrCollection& clusterSums) const {
  // TC distribution
  HGCalTriggerCellSAPtrCollection distributedTCs;
  tcDistribution_.runTriggerCellDistribution(inputs, distributedTCs);

  // Histogramming and seeding
  HGCalHistogramCellSAPtrCollection histogram;
  seeding_.runSeeding(distributedTCs, histogram);

  // Clustering
  HGCalClusterSAPtrCollection protoClusters;
  CentroidHelperPtrCollection readoutFlags;
  clustering_.runClustering(distributedTCs, histogram, clusteredTCs, readoutFlags, protoClusters);

  // Cluster properties
  clusterProperties_.runClusterProperties(protoClusters, readoutFlags, clusterSums);
}
