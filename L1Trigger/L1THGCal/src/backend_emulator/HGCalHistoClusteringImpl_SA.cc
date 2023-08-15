#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringImpl_SA.h"

using namespace std;
using namespace l1thgcfirmware;
HGCalHistoClusteringImplSA::HGCalHistoClusteringImplSA(const ClusterAlgoConfig& config)
    : config_(config), tcDistribution_(config), seeding_(config), clustering_(config), clusterProperties_(config) {}

void HGCalHistoClusteringImplSA::runAlgorithm(const HGCalTriggerCellSAPtrCollections& inputs,
                                              HGCalTriggerCellSAShrPtrCollection& clusteredTCs,
                                              HGCalClusterSAPtrCollection& clusterSums) const {

  // std::cout << "TC Dist" << std::endl;
  // TC distribution
  HGCalTriggerCellSAPtrCollection distributedTCs;
  tcDistribution_.runTriggerCellDistribution(inputs, distributedTCs);

  // std::cout << "Seeding" << std::endl;
  // Histogramming and seeding
  HGCalHistogramCellSAPtrCollection histogram;
  seeding_.runSeeding(distributedTCs, histogram);

  // std::cout << "Clustering" << std::endl;
  // Clustering
  HGCalClusterSAPtrCollection protoClusters;
  CentroidHelperPtrCollection readoutFlags;
  clustering_.runClustering(distributedTCs, histogram, clusteredTCs, readoutFlags, protoClusters);

  // std::cout << "CP" << std::endl;
  // Cluster properties
  clusterProperties_.runClusterProperties(protoClusters, readoutFlags, clusterSums);
  // std::cout << "Done" << std::endl;

}
