#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

class L1GTTRawTrackMVAPatcher : public edm::global::EDProducer<> {
public:
  using L1Track = TTTrack<Ref_Phase2TrackerDigi_>;
  using TTTrackCollection = std::vector<L1Track>;
  using TTTrackCollectionView = edm::View<L1Track>;

  explicit L1GTTRawTrackMVAPatcher(const edm::ParameterSet& iConfig)
      : l1TracksToken_(consumes<TTTrackCollectionView>(iConfig.getParameter<edm::InputTag>("l1TracksInputTag"))),
        outputCollectionName_(iConfig.getParameter<std::string>("outputCollectionName")),
        debug_(iConfig.getParameter<int>("debug")) {
    produces<TTTrackCollection>(outputCollectionName_);
  }

  ~L1GTTRawTrackMVAPatcher() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("l1TracksInputTag", edm::InputTag("l1tTTTracksFromTrackletEmulation", "Level1TTTracks"));
    desc.add<std::string>("outputCollectionName", "Level1TTTracksPatched");
    desc.add<int>("debug", 0);
    descriptions.add("l1GTTRawTrackMVAPatcher", desc);
  }

private:
  void produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup&) const override {
    edm::Handle<TTTrackCollectionView> tracksHandle;
    iEvent.getByToken(l1TracksToken_, tracksHandle);

    auto output = std::make_unique<TTTrackCollection>();
    output->reserve(tracksHandle->size());

    unsigned int itrk = 0;
    for (const auto& inTrack : *tracksHandle) {
      // Start from the original raw track so all existing raw-word fields stay unchanged.
      output->push_back(inTrack);
      auto& patchedTrack = output->back();

      // Temporary copy: let standard CMSSW logic compute the MVA bits from the floating values.
      auto tmpTrack = patchedTrack;
      tmpTrack.setTrackWordBits();

      // Rebuild the final raw word from the ORIGINAL raw fields, but splice in the MVA fields.
      patchedTrack.setTrackWord(
          patchedTrack.getValidWord(),
          patchedTrack.getRinvWord(),
          patchedTrack.getPhiWord(),
          patchedTrack.getTanlWord(),
          patchedTrack.getZ0Word(),
          patchedTrack.getD0Word(),
          patchedTrack.getChi2RPhiWord(),
          patchedTrack.getChi2RZWord(),
          patchedTrack.getBendChi2Word(),
          patchedTrack.getHitPatternWord(),
          tmpTrack.getMVAQualityWord(),
          tmpTrack.getMVAOtherWord());

      if (debug_ >= 2) {
        edm::LogVerbatim("L1GTTRawTrackMVAPatcher")
            << "itrk=" << itrk
            << " rawWordBefore=" << inTrack.getTrackWord().to_string(2)
            << " mvaQ=" << tmpTrack.getMVAQualityWord().to_uint()
            << " mvaOther=" << tmpTrack.getMVAOtherWord().to_uint()
            << " rawWordAfter=" << patchedTrack.getTrackWord().to_string(2);
      }

      ++itrk;
    }

    iEvent.put(std::move(output), outputCollectionName_);
  }

  const edm::EDGetTokenT<TTTrackCollectionView> l1TracksToken_;
  const std::string outputCollectionName_;
  const int debug_;
};

DEFINE_FWK_MODULE(L1GTTRawTrackMVAPatcher);
