#include <algorithm>
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "DataFormats/L1TParticleFlow/interface/PFCluster.h"

#include "L1Trigger/Phase2L1ParticleFlow/interface/deregionizer/deregionizer_input.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/deregionizer/deregionizer_ref.h"

#include "L1Trigger/DemonstratorTools/interface/BoardDataWriter.h"
#include "L1Trigger/DemonstratorTools/interface/utilities.h"

class DeregionizerProducer : public edm::stream::EDProducer<> {
public:
  explicit DeregionizerProducer(const edm::ParameterSet &);
  ~DeregionizerProducer() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  struct LinkWriteInfo {
    l1t::demo::LinkId id;
    size_t payloadWords;
  };

  static std::string interfaceNameForBoard_(uint32_t boardOrder) { return "puppi_in_b" + std::to_string(boardOrder); }

  edm::EDGetTokenT<l1t::PFCandidateRegionalOutput> token_;
  const unsigned int nInputFramesPerBX_;
  l1ct::DeregionizerEmulator emulator_;
  l1ct::DeregionizerInput input_;

  bool writeInputPatternFiles_;
  size_t patternFileBoardTMUX;
  size_t inputGapLength_;

  std::map<std::pair<uint32_t, uint32_t>, LinkWriteInfo> boardLinkToWriteInfo_;
  std::map<l1t::demo::LinkId, size_t> linkPayloadWords_;
  std::map<l1t::demo::LinkId, std::vector<size_t>> channelIdsInput_;
  std::map<std::string, l1t::demo::ChannelSpec> channelSpecsInput_;
  std::unique_ptr<l1t::demo::BoardDataWriter> inputFileWriter_;

  void produce(edm::Event &, const edm::EventSetup &) override;
  void hwToEdm_(const std::vector<l1ct::PuppiObjEmu> &hwOut, std::vector<l1t::PFCandidate> &edmOut) const;
};

DeregionizerProducer::DeregionizerProducer(const edm::ParameterSet &iConfig)
    : token_(consumes<l1t::PFCandidateRegionalOutput>(iConfig.getParameter<edm::InputTag>("RegionalPuppiCands"))),
      nInputFramesPerBX_(iConfig.getParameter<uint32_t>("nInputFramesPerBX")),
      emulator_(iConfig),
      input_(iConfig.getParameter<std::vector<edm::ParameterSet>>("linkConfigs")),
      writeInputPatternFiles_(iConfig.getParameter<bool>("writeInputPatternFiles")),
      patternFileBoardTMUX(0),
      inputGapLength_(0) {
  produces<l1t::PFCandidateCollection>("Puppi");
  produces<l1t::PFCandidateCollection>("TruncatedPuppi");

  if (writeInputPatternFiles_) {
    const auto &pset = iConfig.getParameter<edm::ParameterSet>("inputPatternFilePSet");
    patternFileBoardTMUX = pset.getParameter<uint32_t>("TMUX");
    inputGapLength_ = pset.getParameter<uint32_t>("gapLengthOutput");

    if (patternFileBoardTMUX == 0)
      throw cms::Exception("Configuration") << "inputPatternFilePSet.TMUX must be > 0";

    std::vector<l1ct::DeregionizerInput::BoardInfo> boardInfos = input_.boardInfos_;
    std::sort(boardInfos.begin(), boardInfos.end(), [](const auto &a, const auto &b) { return a.order_ < b.order_; });

    size_t totalLinks = 0;
    size_t maxTmux = patternFileBoardTMUX;
    for (const auto &boardInfo : boardInfos) {
      totalLinks += boardInfo.nLinksPuppi_;
      maxTmux = std::max(maxTmux, static_cast<size_t>(boardInfo.tmuxFactor_));
    }

    if (totalLinks == 0)
      throw cms::Exception("Configuration") << "linkConfigs define zero input links";
    if ((maxTmux % patternFileBoardTMUX) != 0)
      throw cms::Exception("Configuration")
          << "max link tmux " << maxTmux << " is not divisible by inputPatternFilePSet.TMUX=" << patternFileBoardTMUX;

    const size_t nTimeSlices = maxTmux / patternFileBoardTMUX;
    const auto timeSliceConfigs = iConfig.getParameter<std::vector<edm::ParameterSet>>("inputPatternTimeSlices");
    if (timeSliceConfigs.size() != nTimeSlices)
      throw cms::Exception("Configuration")
          << "inputPatternTimeSlices must have " << nTimeSlices << " entries, got " << timeSliceConfigs.size();

    std::vector<std::vector<int32_t>> timeSliceLinks;
    timeSliceLinks.reserve(timeSliceConfigs.size());
    for (size_t iSlice = 0; iSlice < timeSliceConfigs.size(); ++iSlice) {
      const auto links = timeSliceConfigs[iSlice].getParameter<std::vector<int32_t>>("puppiInputLinks");
      if (links.size() != totalLinks)
        throw cms::Exception("Configuration") << "inputPatternTimeSlices[" << iSlice << "].puppiInputLinks has size "
                                              << links.size() << ", expected " << totalLinks;
      timeSliceLinks.push_back(links);
    }

    size_t globalLink = 0;
    for (const auto &boardInfo : boardInfos) {
      if ((boardInfo.tmuxFactor_ % patternFileBoardTMUX) != 0)
        throw cms::Exception("Configuration")
            << "Board order " << boardInfo.order_ << " has tmuxFactor=" << boardInfo.tmuxFactor_
            << " which is not divisible by inputPatternFilePSet.TMUX=" << patternFileBoardTMUX;

      const size_t tmuxRatio = boardInfo.tmuxFactor_ / patternFileBoardTMUX;
      const size_t neededPayloadWords =
          (boardInfo.nPuppiPerRegion_ / boardInfo.nLinksPuppi_) * static_cast<size_t>(boardInfo.regions_.size());
      const size_t maxPayloadWords = boardInfo.tmuxFactor_ * nInputFramesPerBX_;

      if (inputGapLength_ >= maxPayloadWords)
        throw cms::Exception("Configuration")
            << "gapLengthOutput=" << inputGapLength_ << " is too large for board order " << boardInfo.order_
            << " (max payload=" << maxPayloadWords << ")";

      const size_t payloadWords = maxPayloadWords - inputGapLength_;
      if (payloadWords < neededPayloadWords)
        throw cms::Exception("Configuration")
            << "Configured payload words " << payloadWords << " for board order " << boardInfo.order_
            << " is smaller than required words " << neededPayloadWords;

      const std::string interfaceName = interfaceNameForBoard_(boardInfo.order_);
      channelSpecsInput_[interfaceName] = {boardInfo.tmuxFactor_, inputGapLength_, 0};

      for (uint32_t iLink = 0; iLink < boardInfo.nLinksPuppi_; ++iLink, ++globalLink) {
        std::vector<size_t> channelIds;
        channelIds.reserve(tmuxRatio);
        for (size_t iSlice = 0; iSlice < tmuxRatio; ++iSlice) {
          const int32_t channelId = timeSliceLinks[iSlice][globalLink];
          if (channelId < 0)
            throw cms::Exception("Configuration")
                << "Negative channel id " << channelId << " for logical link " << globalLink;
          channelIds.push_back(static_cast<size_t>(channelId));
        }

        l1t::demo::LinkId id{interfaceName, iLink};
        channelIdsInput_[id] = channelIds;
        boardLinkToWriteInfo_[{boardInfo.order_, iLink}] = {id, payloadWords};
        linkPayloadWords_[id] = payloadWords;
      }
    }

    inputFileWriter_ = std::make_unique<l1t::demo::BoardDataWriter>(
        l1t::demo::parseFileFormat(pset.getParameter<std::string>("format")),
        pset.getParameter<std::string>("outputFilename"),
        pset.getParameter<std::string>("outputFileExtension"),
        nInputFramesPerBX_,
        patternFileBoardTMUX,
        pset.getParameter<uint32_t>("maxLinesPerFile"),
        channelIdsInput_,
        channelSpecsInput_);
  }
}

DeregionizerProducer::~DeregionizerProducer() {
  if (inputFileWriter_)
    inputFileWriter_->flush();
}

void DeregionizerProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  auto deregColl = std::make_unique<l1t::PFCandidateCollection>();
  auto truncColl = std::make_unique<l1t::PFCandidateCollection>();

  edm::Handle<l1t::PFCandidateRegionalOutput> src;

  iEvent.getByToken(token_, src);

  std::vector<l1ct::OutputRegion> outputRegions;
  std::vector<l1ct::PuppiObjEmu> hwOut;
  std::vector<l1t::PFCandidate> edmOut;
  std::vector<l1ct::PuppiObjEmu> hwTruncOut;
  std::vector<l1t::PFCandidate> edmTruncOut;

  LogDebug("DeregionizerProducer") << "\nRegional Puppi Candidates";
  for (unsigned int iReg = 0, nReg = src->nRegions(); iReg < nReg; ++iReg) {
    l1ct::OutputRegion tempOutputRegion;

    auto region = src->region(iReg);
    float eta = src->eta(iReg);
    float phi = src->phi(iReg);
    LogDebug("DeregionizerProducer") << "\nRegion " << iReg << "\n"
                                     << "Eta = " << eta << " and Phi = " << phi << "\n"
                                     << "###########";
    for (int i = 0, n = region.size(); i < n; ++i) {
      l1ct::PuppiObjEmu tempPuppi;
      const l1t::PFCandidate &cand = region[i];

      tempPuppi.initFromBits(cand.encodedPuppi64());
      tempPuppi.srcCand = &cand;
      tempOutputRegion.puppi.push_back(tempPuppi);
      LogDebug("DeregionizerProducer") << "pt[" << i << "] = " << tempOutputRegion.puppi.back().hwPt << ", eta[" << i
                                       << "] = " << tempOutputRegion.puppi.back().floatEta() << ", phi[" << i
                                       << "] = " << tempOutputRegion.puppi.back().floatPhi();
    }
    outputRegions.push_back(tempOutputRegion);
  }

  const auto layer2InWithPlacement = input_.orderInputsWithPlacement(outputRegions);
  std::vector<std::vector<std::vector<l1ct::PuppiObjEmu>>> layer2In(layer2InWithPlacement.size());
  for (size_t iClock = 0; iClock < layer2InWithPlacement.size(); ++iClock) {
    layer2In[iClock].resize(layer2InWithPlacement[iClock].size());
    for (size_t iBoard = 0; iBoard < layer2InWithPlacement[iClock].size(); ++iBoard) {
      for (const auto &placedPuppi : layer2InWithPlacement[iClock][iBoard]) {
        layer2In[iClock][iBoard].push_back(placedPuppi.first);
      }
    }
  }

  if (writeInputPatternFiles_) {
    std::map<l1t::demo::LinkId, std::vector<ap_uint<64>>> links;
    for (const auto &[id, payloadWords] : linkPayloadWords_)
      links.emplace(id, std::vector<ap_uint<64>>(payloadWords, ap_uint<64>(0)));

    for (const auto &clockSlice : layer2InWithPlacement) {
      for (const auto &boardSlice : clockSlice) {
        for (const auto &entry : boardSlice) {
          const auto &obj = entry.first;
          const auto &lpi = entry.second;
          auto it = boardLinkToWriteInfo_.find({lpi.board_, lpi.link_});
          if (it == boardLinkToWriteInfo_.end())
            continue;
          const auto &info = it->second;
          if (lpi.clock_cycle_ < info.payloadWords) {
            links[info.id][lpi.clock_cycle_] = obj.pack();
          }
        }
      }
    }
    l1t::demo::EventData eventDataInputs;
    for (const auto &[id, words] : links)
      eventDataInputs.add(id, words);
    inputFileWriter_->addEvent(eventDataInputs);
  }

  emulator_.run(layer2In, hwOut, hwTruncOut);

  DeregionizerProducer::hwToEdm_(hwOut, edmOut);
  DeregionizerProducer::hwToEdm_(hwTruncOut, edmTruncOut);

  deregColl->swap(edmOut);
  truncColl->swap(edmTruncOut);

  iEvent.put(std::move(deregColl), "Puppi");
  iEvent.put(std::move(truncColl), "TruncatedPuppi");
}

void DeregionizerProducer::hwToEdm_(const std::vector<l1ct::PuppiObjEmu> &hwOut,
                                    std::vector<l1t::PFCandidate> &edmOut) const {
  for (const auto &hwPuppi : hwOut) {
    l1t::PFCandidate::ParticleType type;
    float mass = 0.13f;
    if (hwPuppi.hwId.charged()) {
      if (hwPuppi.hwId.isMuon()) {
        type = l1t::PFCandidate::Muon;
        mass = 0.105;
      } else if (hwPuppi.hwId.isElectron()) {
        type = l1t::PFCandidate::Electron;
        mass = 0.005;
      } else
        type = l1t::PFCandidate::ChargedHadron;
    } else {
      type = hwPuppi.hwId.isPhoton() ? l1t::PFCandidate::Photon : l1t::PFCandidate::NeutralHadron;
      mass = hwPuppi.hwId.isPhoton() ? 0.0 : 0.5;
    }
    reco::Particle::PolarLorentzVector p4(hwPuppi.floatPt(), hwPuppi.floatEta(), hwPuppi.floatPhi(), mass);
    edmOut.emplace_back(
        type, hwPuppi.intCharge(), p4, hwPuppi.floatPuppiW(), hwPuppi.intPt(), hwPuppi.intEta(), hwPuppi.intPhi());
    if (hwPuppi.hwId.charged()) {
      edmOut.back().setZ0(hwPuppi.floatZ0());
      edmOut.back().setDxy(hwPuppi.floatDxy());
      edmOut.back().setHwZ0(hwPuppi.hwZ0());
      edmOut.back().setHwDxy(hwPuppi.hwDxy());
      edmOut.back().setHwTkQuality(hwPuppi.hwTkQuality());
    } else {
      edmOut.back().setHwPuppiWeight(hwPuppi.hwPuppiW());
      edmOut.back().setHwEmID(hwPuppi.hwEmID());
    }
    edmOut.back().setEncodedPuppi64(hwPuppi.pack().to_uint64());
    edmOut.back().setCaloPtr(hwPuppi.srcCand->caloPtr());
    edmOut.back().setPFTrack(hwPuppi.srcCand->pfTrack());
    edmOut.back().setMuon(hwPuppi.srcCand->muon());
  }
}

void DeregionizerProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("RegionalPuppiCands", edm::InputTag("l1tLayer1", "PuppiRegional"));
  desc.add<unsigned int>("nPuppiFinalBuffer", 128);
  desc.add<unsigned int>("nPuppiPerClk", 6);
  desc.add<unsigned int>("nPuppiFirstBuffers", 12);
  desc.add<unsigned int>("nPuppiSecondBuffers", 32);
  desc.add<unsigned int>("nPuppiThirdBuffers", 64);
  desc.add<unsigned int>("nInputFramesPerBX", 9);
  desc.add<bool>("writeInputPatternFiles", false);

  edm::ParameterSetDescription inputPatternPSet;
  inputPatternPSet.add<uint32_t>("gapLengthOutput", 0);
  inputPatternPSet.add<uint32_t>("TMUX", 6);
  inputPatternPSet.add<uint32_t>("maxLinesPerFile", 1024);
  inputPatternPSet.add<std::string>("outputFilename", "L1DeregionizerInput");
  inputPatternPSet.add<std::string>("format", "EMPv2");
  inputPatternPSet.add<std::string>("outputFileExtension", "txt.gz");
  desc.add<edm::ParameterSetDescription>("inputPatternFilePSet", inputPatternPSet);

  edm::ParameterSetDescription inputTimeSlicePSet;
  inputTimeSlicePSet.add<std::vector<int32_t>>("puppiInputLinks", {});
  desc.addVPSet("inputPatternTimeSlices", inputTimeSlicePSet, std::vector<edm::ParameterSet>{});

  edm::ParameterSetDescription linkConfigDummyValidator;
  linkConfigDummyValidator.setAllowAnything();
  desc.addVPSet("linkConfigs", linkConfigDummyValidator);
  descriptions.add("DeregionizerProducer", desc);
}

DEFINE_FWK_MODULE(DeregionizerProducer);
