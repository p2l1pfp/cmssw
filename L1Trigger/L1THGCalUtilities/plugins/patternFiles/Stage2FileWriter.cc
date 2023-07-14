// -*- C++ -*-
//
// Package:    L1Trigger/L1THGCalUtilities
// Class:      Stage2FileWriter
//
/**\class Stage2FileWriter Stage2FileWriter.cc L1Trigger/L1THGCalUtilities/plugins/patternFiles/Stage2FileWriter.cc

 Description: EDAnalyzer for writing I/O buffer files for hardware/firmware tests of HGC stage 2

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
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"
#include "DataFormats/L1THGCal/interface/HGCalCluster_HW.h"
#include "DataFormats/Common/interface/View.h"

#include "L1Trigger/DemonstratorTools/interface/BoardDataWriter.h"
#include "L1Trigger/DemonstratorTools/interface/utilities.h"

#include "ap_int.h"

//
// class declaration
//

class Stage2FileWriter : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit Stage2FileWriter(const edm::ParameterSet&);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // ----------constants, enums and typedefs ---------
  static constexpr size_t kFramesPerTMUXPeriod = 9;
  static constexpr size_t kGapLengthOutput = 0;
  static constexpr size_t kS2BoardTMUX = 18;
  static constexpr size_t kMaxLinesPerFile = 648;

  const std::map<l1t::demo::LinkId, std::pair<l1t::demo::ChannelSpec, std::vector<size_t>>>
      kChannelSpecsOutputToL1T = {
          /* logical channel within time slice -> {{link TMUX, inter-packet gap}, vector of channel indices} */
          {{"towersAndClusters", 0}, {{kS2BoardTMUX, kGapLengthOutput}, {81}}},
          {{"towersAndClusters", 1}, {{kS2BoardTMUX, kGapLengthOutput}, {82}}},
          {{"towersAndClusters", 2}, {{kS2BoardTMUX, kGapLengthOutput}, {83}}},
          {{"towersAndClusters", 3}, {{kS2BoardTMUX, kGapLengthOutput}, {84}}}
        };

  const std::map<l1t::demo::LinkId, std::pair<l1t::demo::ChannelSpec, std::vector<size_t>>>
      kChannelSpecsClusterSumsInput = {
          /* logical channel within time slice -> {{link TMUX, inter-packet gap}, vector of channel indices} */
          {{"clusterSumRecord", 0}, {{kS2BoardTMUX, kGapLengthOutput}, {80}}},
          {{"clusterSumRecord", 1}, {{kS2BoardTMUX, kGapLengthOutput}, {81}}},
          {{"clusterSumRecord", 2}, {{kS2BoardTMUX, kGapLengthOutput}, {82}}},
          {{"clusterSumRecord", 3}, {{kS2BoardTMUX, kGapLengthOutput}, {83}}},
          {{"clusterSumRecord", 4}, {{kS2BoardTMUX, kGapLengthOutput}, {84}}},
          {{"clusterSumRecord", 5}, {{kS2BoardTMUX, kGapLengthOutput}, {85}}},
          {{"clusterSumRecord", 6}, {{kS2BoardTMUX, kGapLengthOutput}, {86}}},
          {{"clusterSumRecord", 7}, {{kS2BoardTMUX, kGapLengthOutput}, {87}}},
        };

  // typedef TTTrack<Ref_Phase2TrackerDigi_> Track_t;

  // ----------member functions ----------------------
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  std::array<std::vector<ap_uint<64>>, 4> encodeTowersAndClusters(const l1t::HGCalMulticlusterBxCollection&, const unsigned int iSector );

  std::array<std::vector<ap_uint<64>>, 8> encodeClusterSumRecord(const l1t::HGCalMulticlusterBxCollection&, const unsigned int iSector );

  // ----------member data ---------------------------
  edm::EDGetTokenT<l1t::HGCalMulticlusterBxCollection> clustersToken_;

  std::vector<l1t::demo::BoardDataWriter> fileWritersOutputToL1T_;
  std::vector<l1t::demo::BoardDataWriter> fileWritersClusterSumsInput_;
  unsigned int tmIndex_;
  unsigned int eventCounter_;
};

//
// class implementation
//

Stage2FileWriter::Stage2FileWriter(const edm::ParameterSet& iConfig)
    : clustersToken_(consumes<l1t::HGCalMulticlusterBxCollection>(iConfig.getUntrackedParameter<edm::InputTag>("clusters"))),
      tmIndex_(iConfig.getUntrackedParameter<unsigned int>("tmIndex")),
      eventCounter_(0) {
      for ( unsigned int iFileWriter=0; iFileWriter < 6; ++iFileWriter ) {
        fileWritersOutputToL1T_.emplace_back(l1t::demo::parseFileFormat(iConfig.getUntrackedParameter<std::string>("format")),
                                    std::string("hgc_sec")+std::to_string(iFileWriter)+"_tm"+std::to_string(tmIndex_)+std::string("-output-ref"),
                                    kFramesPerTMUXPeriod,
                                    kS2BoardTMUX,
                                    kMaxLinesPerFile,
                                    kChannelSpecsOutputToL1T);
        fileWritersClusterSumsInput_.emplace_back(l1t::demo::parseFileFormat(iConfig.getUntrackedParameter<std::string>("format")),
                                            std::string("hgc_sec")+std::to_string(iFileWriter)+"_tm"+std::to_string(tmIndex_)+std::string("-input"),
                                            kFramesPerTMUXPeriod,
                                            kS2BoardTMUX,
                                            kMaxLinesPerFile,
                                            kChannelSpecsClusterSumsInput);
      }
    }

void Stage2FileWriter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  // Loop S2 sectors and produce pattern file for each sector
  // Would be better to sort clusters into sectors first
  for ( unsigned int iSector = 0; iSector < 6; ++iSector ) {
    // 1) Encode tower and cluster information onto vectors containing link data
    const auto outputData(encodeTowersAndClusters(iEvent.get(clustersToken_), iSector) );
    const auto clusterSumInputData(encodeClusterSumRecord(iEvent.get(clustersToken_), iSector) );

    // 2) Pack track information into 'event data' object, and pass that to file writer
    l1t::demo::EventData eventDataTowersAndClusters;
    for (size_t i = 0; i < 4; i++) {
      eventDataTowersAndClusters.add({"towersAndClusters", i}, outputData.at(i));
    }
    fileWritersOutputToL1T_.at(iSector).addEvent(eventDataTowersAndClusters);

    l1t::demo::EventData eventDataClusterSums;
    for (size_t i = 0; i < 8; i++) {
      eventDataClusterSums.add({"clusterSumRecord", i}, clusterSumInputData.at(i));
    }
    fileWritersClusterSumsInput_.at(iSector).addEvent(eventDataClusterSums);
  }
  ++eventCounter_;
}

// ------------ method called once each job just after ending the event loop  ------------
void Stage2FileWriter::endJob() {
  // Writing pending events to file before exiting
  for ( unsigned int iSector = 0; iSector < 6; ++iSector ) {
    fileWritersOutputToL1T_.at(iSector).flush();
    fileWritersClusterSumsInput_.at(iSector).flush();
  }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Stage2FileWriter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // Stage2FileWriter
  edm::ParameterSetDescription desc;
  desc.addUntracked<edm::InputTag>("clusters", edm::InputTag("hgcalBackEndLayer2Producer", "HGCalBackendLayer2Processor3DClusteringSA"));
  desc.addUntracked<unsigned int>("tmIndex", 0);
  desc.addUntracked<std::string>("format", "EMP");
  descriptions.add("Stage2FileWriter", desc);
}

std::array<std::vector<ap_uint<64>>, 4> Stage2FileWriter::encodeTowersAndClusters(const l1t::HGCalMulticlusterBxCollection& clusters, const unsigned int iSector ) {
  std::array<std::vector<ap_uint<64>>, 4> output;

  // Decode which zside and sector we are extracting data for
  int zside = (iSector > 2) ? 1 : -1;
  unsigned int sector = iSector % 3;
  if ( zside == 1 ) {
    if ( sector == 1 ) sector = 2;
    else if ( sector == 2 ) sector = 1;
  }
  // std::cout << "Sectors : " << iSector << " " << zside << " " << sector << std::endl;

  // First frame empty for alignment
  ap_uint<64> packetHeader = 0;
  packetHeader(63,48) = 0xabc0;
  packetHeader(31,24) = iSector;
  for ( unsigned iLink = 0; iLink<4; ++iLink ) {
    packetHeader(23,16) = iLink;
    output[iLink].push_back( packetHeader );
  }

  // Tower data for this sector
  // Dummy data
  for ( unsigned int iTowerFrame = 0; iTowerFrame < 30; ++iTowerFrame ) {
    output[0].push_back( iTowerFrame*4 + 0 );
    output[1].push_back( iTowerFrame*4 + 1 );
    output[2].push_back( iTowerFrame*4 + 2 );
    output[3].push_back( iTowerFrame*4 + 3 );
  }

  // Cluster data for this sector
  unsigned int iCluster = 0;
  for (auto cl3d_itr = clusters.begin(0); cl3d_itr != clusters.end(0); cl3d_itr++) {
    if ( cl3d_itr->getHwZSide() != zside ) continue;
    if ( cl3d_itr->getHwSector() != sector ) continue;

    ++iCluster;
    if ( iCluster > 160 ) break;
    // if ( sector == 0 && zside == -1 ) {
    //   std::cout << "Cluster pt, eta, phi, nTC,  : " << sector << " " << cl3d_itr->pt() << " " << cl3d_itr->eta() << " " << cl3d_itr->phi() << " " << cl3d_itr->constituents().size() << " " << cl3d_itr->hw_sigma_e_quotient() << " " << cl3d_itr->hw_sigma_e_fraction() << std::endl;
    // }
    //  << " " << cl3d_itr->getHwData()[0].to_string() << std::endl;
    // std::cout << "Phi, eta : " << cl3d_itr->phi() << " " << cl3d_itr->eta() << " " << cl3d_itr->getHwData()[1].to_string() << std::endl;
    const auto& clusterWords = cl3d_itr->getHwData();
    // std::cout << "Cluster words : " ;
    // for (const auto& word : clusterWords ) std::cout << word << " ";
    // std::cout << std::endl;
    output[0].push_back( clusterWords[0] );
    output[1].push_back( clusterWords[1] );
    output[2].push_back( clusterWords[2] );
    output[3].push_back( clusterWords[3] );

    // if ( sector == 0 && zside == -1 ) {

    //   const auto unpacked = l1thgcfirmware::HGCalCluster_HW::unpack(cl3d_itr->getHwData());
    //   if ( abs( cl3d_itr->pt() - unpacked.e*0.25 ) > 0.001 ||
    //       abs( cl3d_itr->iPt(l1t::HGCalMulticluster::EnergyInterpretation::EM) - unpacked.e_em*0.25 ) > 0.001 ||
    //       abs( abs(cl3d_itr->eta()) - unpacked.w_eta * M_PI/720 ) > 0.001 ||
    //       abs( cl3d_itr->phi() - ( unpacked.w_phi * M_PI/720 + M_PI / 2 ) ) > 0.001 ||
    //       abs( cl3d_itr->zBarycenter() - unpacked.w_z * 0.05 ) > 0.001 ||
    //       abs( cl3d_itr->sigmaRRTot() - unpacked.sigma_roz * 0.0001920625 ) > 0.001
    //   ) {
    //     std::cout << "---> PANIC : Comparing unpacked hw data with edm object" << std::endl;
    //     std::cout << cl3d_itr->pt() << " " << unpacked.e*0.25 << std::endl;
    //     std::cout << cl3d_itr->iPt(l1t::HGCalMulticluster::EnergyInterpretation::EM) << " " << unpacked.e_em*0.25 << std::endl;
    //     std::cout << abs(cl3d_itr->eta()) << " " << unpacked.w_eta * M_PI/720 << std::endl;
    //     std::cout << cl3d_itr->phi() << " " << unpacked.w_phi * M_PI/720 + M_PI / 2 << std::endl;
    //     std::cout << cl3d_itr->zBarycenter() << " " << unpacked.w_z * 0.05 << std::endl;
    //     std::cout << cl3d_itr->sigmaRRTot() << " " << unpacked.sigma_roz * 0.0001920625 << std::endl;
    //   }
    // }
    
  }

  return output;
}

std::array<std::vector<ap_uint<64>>, 8> Stage2FileWriter::encodeClusterSumRecord(const l1t::HGCalMulticlusterBxCollection& clusters, const unsigned int iSector ) {

  std::array<std::vector<ap_uint<64>>, 8> output;

  // Decode which zside and sector we are extracting data for
  int zside = (iSector > 2) ? 1 : -1;
  unsigned int sector = iSector % 3;
  if ( zside == 1 ) {
    if ( sector == 1 ) sector = 2;
    else if ( sector == 2 ) sector = 1;
  }
  // std::cout << "Sectors : " << iSector << " " << zside << " " << sector << std::endl;

  ap_uint<16> eventCount = eventCounter_;
  ap_uint<8> boardID = iSector;
  ap_uint<8> tmIndex = 0;

  ap_uint<64> headerWord = (eventCount, boardID, tmIndex, ap_uint<32>(0));

  unsigned int nClusterSums = 0;

  for (auto cl3d_itr = clusters.begin(0); cl3d_itr != clusters.end(0); cl3d_itr++) {
    if ( cl3d_itr->getHwZSide() != zside ) continue;
    if ( cl3d_itr->getHwSector() != sector ) continue;

    // ++iCluster;
    // if ( iCluster > 160 ) break;
    // std::cout << "Adding cluster, sector : " << iSector << " " << zside << " " << sector << " " << cl3d_itr->pt() << " " << cl3d_itr->eta() << " " << cl3d_itr->phi() << " " << cl3d_itr->size() << std::endl;
    //  << " " << cl3d_itr->getHwData()[0].to_string() << std::endl;
    // std::cout << "Phi, eta : " << cl3d_itr->phi() << " " << cl3d_itr->eta() << " " << cl3d_itr->getHwData()[1].to_string() << std::endl;
    const auto& clusterSumWords = cl3d_itr->getHwClusterSumData();
    for ( unsigned iWord = 0; iWord < 7; ++iWord ) {
      // std::cout << "Cluster sum word : " << iWord << " " << clusterSumWords[iWord] << std::endl;
      output[iWord].push_back( clusterSumWords[iWord] );
    }
    output[7].push_back( headerWord );
    ++nClusterSums;
  }
  // std::cout << "Added all clusters" << std::endl;

  // Add dummy entry if there weren't any cluster sums in this sector in this event
  if ( output[0].size() == 0 ) {
    for ( unsigned iWord = 0; iWord < 7; ++iWord ) {
      output[iWord].push_back( 0 );
    }
    output[7].push_back( headerWord );
    ++nClusterSums;
  }

  while ( nClusterSums < 161 ) {
    // for ( unsigned iWord = 0; iWord < 7; ++iWord ) {
    //   output[iWord].push_back( 0 );
    // }
    output[7].push_back( headerWord );
    ++nClusterSums;
  }
  return output;
}


//define this as a plug-in
DEFINE_FWK_MODULE(Stage2FileWriter);
