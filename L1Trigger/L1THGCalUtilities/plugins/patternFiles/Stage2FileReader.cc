// -*- C++ -*-
//
// Package:    L1Trigger/L1THGCalUtilities
// Class:      Stage2FileReader
//
/**\class Stage2FileReader Stage2FileReader.cc L1Trigger/L1THGCalUtilities/plugins/patternFiles/Stage2FileReader.cc

 Description: EDAnalyzer for reading I/O buffer files for hardware/firmware tests of HGC stage 2

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Emyr Clement
//         Created:  Tue, 13 Apr 2022
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"
#include "DataFormats/L1THGCal/interface/HGCalCluster_HW.h"
#include "L1Trigger/DemonstratorTools/interface/BoardDataReader.h"
#include "L1Trigger/DemonstratorTools/interface/utilities.h"
#include "L1Trigger/L1THGCalUtilities/interface/patternFiles/codecs_clusters.h"

//
// class declaration
//

class Stage2FileReader : public edm::stream::EDProducer<> {
public:
  explicit Stage2FileReader(const edm::ParameterSet&);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // ----------constants, enums and typedefs ---------
  static constexpr size_t kFramesPerTMUXPeriod = 9;
  static constexpr size_t kGapLengthOutput = 0;
  static constexpr size_t kS2BoardTMUX = 18;
  static constexpr size_t kEmptyFrames = 0;

  const std::map<l1t::demo::LinkId, std::pair<l1t::demo::ChannelSpec, std::vector<size_t>>>
      kChannelSpecsOutputToL1T = {
          /* logical channel within time slice -> {{link TMUX, inter-packet gap}, vector of channel indices} */
          {{"towersAndClusters", 0}, {{kS2BoardTMUX, kGapLengthOutput}, {0}}},
          {{"towersAndClusters", 1}, {{kS2BoardTMUX, kGapLengthOutput}, {1}}},
          {{"towersAndClusters", 2}, {{kS2BoardTMUX, kGapLengthOutput}, {2}}},
          {{"towersAndClusters", 3}, {{kS2BoardTMUX, kGapLengthOutput}, {3}}}
        };

  // ----------member functions ----------------------
  void produce(edm::Event&, const edm::EventSetup&) override;
  std::vector<l1thgcfirmware::HGCalCluster_HW> getRefHWClusters( const l1t::HGCalMulticlusterBxCollection& refClusters, const unsigned int iSector );
  void compareClustersToRef( std::vector<l1thgcfirmware::HGCalCluster_HW> clusters, std::vector<l1thgcfirmware::HGCalCluster_HW> refClusters );

  // ----------member data ---------------------------
  l1t::demo::BoardDataReader fileReader_;

  edm::EDGetTokenT<l1t::HGCalMulticlusterBxCollection> refClustersToken_;

  unsigned int nClusters_;

};

//
// class implementation
//

Stage2FileReader::Stage2FileReader(const edm::ParameterSet& iConfig)
    : fileReader_(l1t::demo::parseFileFormat(iConfig.getUntrackedParameter<std::string>("format")),
                  iConfig.getParameter<std::vector<std::string>>("files"),
                  kFramesPerTMUXPeriod,
                  kS2BoardTMUX,
                  kEmptyFrames,
                  kChannelSpecsOutputToL1T),
    refClustersToken_(consumes<l1t::HGCalMulticlusterBxCollection>(iConfig.getUntrackedParameter<edm::InputTag>("refClustersTag"))),
    nClusters_(0)
 {
  produces<l1t::HGCalMulticlusterBxCollection>();
}

// ------------ method called to produce the data  ------------
void Stage2FileReader::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
//   using namespace l1t::demo::codecs;

  l1t::demo::EventData eventData(fileReader_.getNextEvent());

  const unsigned bitsPerWord = 64;
  const unsigned wordsPerCluster = 4;
  std::array<std::vector<ap_uint<bitsPerWord>>, wordsPerCluster> clusterWords;
  for ( unsigned iWord = 0; iWord < wordsPerCluster; ++iWord ) {
      clusterWords.at(iWord) = eventData.at({"towersAndClusters", iWord});
  }


  std::vector<l1thgcfirmware::HGCalCluster_HW> hwClusters = l1thgcfirmware::decodeClusters(clusterWords);

  // Ref clusters
  const l1t::HGCalMulticlusterBxCollection refClusters = iEvent.get(refClustersToken_);
  std::vector<l1thgcfirmware::HGCalCluster_HW> refHWClusters = getRefHWClusters( refClusters, 0 );

  compareClustersToRef( hwClusters, refHWClusters );
  std::cout << "N clusters so far : " << nClusters_ << std::endl;
}

std::vector<l1thgcfirmware::HGCalCluster_HW> Stage2FileReader::getRefHWClusters( const l1t::HGCalMulticlusterBxCollection& refClusters, const unsigned int iSector ) {

  // Decode which zside and sector we are extracting data for
  int zside = (iSector > 2) ? 1 : -1;
  unsigned int sector = iSector % 3;
  if ( zside == 1 ) {
    if ( sector == 1 ) sector = 2;
    else if ( sector == 2 ) sector = 1;
  }

  std::vector<l1thgcfirmware::HGCalCluster_HW> refHWClusters;
  unsigned int iCluster = 0;
  for (auto cl3d_itr = refClusters.begin(0); cl3d_itr != refClusters.end(0); cl3d_itr++) {
    if ( cl3d_itr->getHwZSide() != zside ) continue;
    if ( cl3d_itr->getHwSector() != sector ) continue;

    ++iCluster;
    if ( iCluster > 160 ) break;

    const auto& clusterWords = cl3d_itr->getHwData();
    refHWClusters.emplace_back( l1thgcfirmware::HGCalCluster_HW::unpack(clusterWords) );
  }
  return refHWClusters;
}

void Stage2FileReader::compareClustersToRef( std::vector<l1thgcfirmware::HGCalCluster_HW> clusters, std::vector<l1thgcfirmware::HGCalCluster_HW> refClusters ) {
  nClusters_ += refClusters.size();
  if ( clusters.size() != refClusters.size() ) {
    std::cout << "---> Different number of clusters : " << refClusters.size() << " " << clusters.size() << std::endl;
    // return;
  }

  bool allGood = true;
  for ( unsigned int iCluster = 0; iCluster < refClusters.size(); ++iCluster ) { // BUG IN FIRMWARE!!!  First event has duplicate cluster, all others missing last cluster...
    if ( iCluster >= clusters.size() ) {
      std::cout << "Skipping extra ref cluster : " << iCluster << std::endl;
      break;
    }
    const auto& refCluster = refClusters.at(iCluster);
    const auto& cluster = clusters.at(iCluster);  // BUG IN FIRMWARE!!!  First event has duplicate cluster, all others missing last cluster...
    if ( cluster != refCluster ) {
      std::cout << "Cluster comparison : " << iCluster << std::endl;
      // Catch common/known issues
      if ( cluster.sigma_roz == 127 && refCluster.sigma_roz != 127 ) {
        std::cout << "Saturated sigma r/z r/z in firmware" << std::endl;
        allGood = false;
        continue;
      }
      if ( abs(refCluster.w_eta - cluster.w_eta) == 1 ) {
        std::cout << "Etas differ by one bit" << std::endl;
        allGood = false;
        continue;
      }
      if ( clusters.size() > iCluster+1) {
        if ( refCluster == clusters.at(iCluster+1 )) {
          std::cout << "Cluster ordering out by one" << std::endl;
          continue;
        }
      }
      std::cout << "--->  e : " << refCluster.e << " " << cluster.e << std::endl;
      std::cout << "--->  e EM : " << refCluster.e_em << " " << cluster.e_em << std::endl;
      std::cout << "--->  eta : " << refCluster.w_eta << " " << cluster.w_eta << std::endl;
      std::cout << "--->  phi : " << refCluster.w_phi << " " << cluster.w_phi << std::endl;
      std::cout << "--->  z : " << refCluster.w_z << " " << cluster.w_z << std::endl;
      std::cout << "--->  sigma r/z : " << refCluster.sigma_roz << " " << cluster.sigma_roz << std::endl;
      std::cout << "-----------------------------------" << std::endl;
      allGood = false;
    }
  }

  if ( allGood ) std::cout << "---> Good event, all match" << std::endl;

  // for (const auto& refCluster : refClusters ) {
  //   std::cout << "Ref cluster : " << refCluster.e << " " << refCluster.w_z << std::endl;
  // }
  // for (const auto& cluster : clusters ) {
  //   if (cluster.e == 0 ) continue;
  //   std::cout << "Cluster : " << cluster.e << " " << cluster.w_z << std::endl;
  // }

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Stage2FileReader::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // Stage2FileReader
  edm::ParameterSetDescription desc;
  desc.add<std::vector<std::string>>("files",
                                     {
                                         "HGCS2OutputToL1TFile_Sector0_0.txt",
                                     });
  desc.addUntracked<std::string>("format", "EMPv2");
  desc.addUntracked<edm::InputTag>("refClustersTag", edm::InputTag("l1tHGCalBackEndLayer2Producer", "HGCalBackendLayer2Processor3DClusteringSA"));
  descriptions.add("Stage2FileReader", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Stage2FileReader);