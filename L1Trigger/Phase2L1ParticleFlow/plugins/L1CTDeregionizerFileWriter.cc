#include <memory>
#include <numeric>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/View.h"

#include "L1Trigger/DemonstratorTools/interface/BoardDataWriter.h"
#include "L1Trigger/DemonstratorTools/interface/utilities.h"
#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"

//
// class declaration
//

class L1CTDeregionizerFileWriter : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit L1CTDeregionizerFileWriter(const edm::ParameterSet&);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // ----------constants, enums and typedefs ---------
  edm::EDGetTokenT<std::vector<l1t::PFCandidate>> l1PFToken_;
  size_t nFramesPerBX_;
  size_t ctl2BoardTMUX_;
  size_t gapLengthOutput_;
  size_t maxLinesPerFile_;
  std::map<l1t::demo::LinkId, std::pair<l1t::demo::ChannelSpec, std::vector<size_t>>> channelSpecsOutput_;

  // ----------member functions ----------------------
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  std::vector<std::vector<ap_uint<64>>> encodePuppi(const std::vector<l1t::PFCandidate> particles, const int nrows);

  l1t::demo::BoardDataWriter fileWriterOutput_;
};

L1CTDeregionizerFileWriter::L1CTDeregionizerFileWriter(const edm::ParameterSet& iConfig)
    : l1PFToken_(consumes<std::vector<l1t::PFCandidate>>(iConfig.getParameter<edm::InputTag>("L1PFObjects"))),
      nFramesPerBX_(iConfig.getParameter<unsigned>("nFramesPerBX")),
      ctl2BoardTMUX_(iConfig.getParameter<unsigned>("TMUX")),
      gapLengthOutput_(iConfig.getParameter<unsigned>("gapLengthOutput")),
      maxLinesPerFile_(iConfig.getParameter<unsigned>("maxLinesPerFile")),
      channelSpecsOutput_{
          {{"deregionizer", 0}, {{ctl2BoardTMUX_, gapLengthOutput_}, {0}}},
          {{"deregionizer", 1}, {{ctl2BoardTMUX_, gapLengthOutput_}, {1}}},
          {{"deregionizer", 2}, {{ctl2BoardTMUX_, gapLengthOutput_}, {2}}},
          {{"deregionizer", 3}, {{ctl2BoardTMUX_, gapLengthOutput_}, {3}}},
      },
      fileWriterOutput_(l1t::demo::parseFileFormat(iConfig.getParameter<std::string>("format")),
                        iConfig.getParameter<std::string>("outputFilename"),
                        iConfig.getParameter<std::string>("outputFileExtension"),
                        nFramesPerBX_,
                        ctl2BoardTMUX_,
                        maxLinesPerFile_,
                        channelSpecsOutput_) {}

void L1CTDeregionizerFileWriter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  const std::vector<l1t::PFCandidate> particles = iEvent.get(l1PFToken_);
  std::vector<std::vector<ap_uint<64>>> link_words = encodePuppi(particles, 4);

  l1t::demo::EventData eventData;
  for (size_t i = 0; i < 4; i++) {
    eventData.add({"deregionizer", i}, link_words.at(i));
  }
  fileWriterOutput_.addEvent(eventData);
}

// ------------ method called once each job just after ending the event loop  ------------
void L1CTDeregionizerFileWriter::endJob() {
  // Writing pending events to file before exiting
  fileWriterOutput_.flush();
}

std::vector<std::vector<ap_uint<64>>> L1CTDeregionizerFileWriter::encodePuppi(
    const std::vector<l1t::PFCandidate> particles, const int nrows) {
  // 'reshape' the 1D View of particles to a 2D vector with nrows in the first dimension
  // pack the particles to their 64 bit HW representation
  std::vector<std::vector<ap_uint<64>>> particles_packed_reshaped(nrows);
  for (uint i = 0; i < particles.size(); i++) {
    ap_uint<64> p = particles.at(i).encodedPuppi64();
    particles_packed_reshaped[i % nrows].push_back(p);
  }
  return particles_packed_reshaped;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void L1CTDeregionizerFileWriter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("L1PFObjects", edm::InputTag("l1tLayer2Deregionizer", "Puppi"));
  desc.add<std::string>("outputFilename", "L1CTDeregionizerPatterns");
  desc.add<std::string>("outputFileExtension", "txt.gz");
  desc.add<uint32_t>("nFramesPerBX", 9);
  desc.add<uint32_t>("gapLengthOutput", 4);
  desc.add<uint32_t>("TMUX", 6);
  desc.add<uint32_t>("maxLinesPerFile", 1024);
  desc.add<std::string>("format", "EMPv2");
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1CTDeregionizerFileWriter);
