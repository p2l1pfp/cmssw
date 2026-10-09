#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/FileInPath.h"

#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/L1TEmbeddingModelEmulator.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/L1TEmbeddingModelServer.h"

#include <iostream>
#include <memory>
#include <vector>

class L1TEmbeddingModelProducer : public edm::global::EDProducer<> {
public:
  explicit L1TEmbeddingModelProducer(const edm::ParameterSet &);
  ~L1TEmbeddingModelProducer() override = default;
  static void fillDescriptions(edm::ConfigurationDescriptions &);

private:
  void produce(edm::StreamID, edm::Event &, const edm::EventSetup &) const override;

  // Model name passed to hls4mlEmulator::ModelLoader, which dlopens <name>.so.
  //   Local test: package-relative path WITH "/", e.g. "L1Trigger/Phase2L1ParticleFlow/data/ctl2_embedding_v0"
  //               -> resolved via edm::FileInPath so it works regardless of cwd.
  //   Production: bare name WITHOUT "/", e.g. "L1TEmbeddingModel"
  //               -> found on LD_LIBRARY_PATH from the L1TEmbeddingModel RPM external.
  static std::string resolveModelName(const std::string &name);

  edm::EDGetTokenT<std::vector<l1t::PFCandidate>> puppiToken_;
  bool saveInput_;
  std::unique_ptr<L1TEmbeddingModelServer> model_;
};

std::string L1TEmbeddingModelProducer::resolveModelName(const std::string &name) {
  if (name.empty() || name.find('/') == std::string::npos)
    return name;  // bare name: resolved from LD_LIBRARY_PATH in production
  const std::string so = name + ".so";
  const std::string full = edm::FileInPath(so).fullPath();
  return full.substr(0, full.size() - 3);
}

L1TEmbeddingModelProducer::L1TEmbeddingModelProducer(const edm::ParameterSet &cfg)
    : puppiToken_(consumes<std::vector<l1t::PFCandidate>>(cfg.getParameter<edm::InputTag>("PuppiCandidates"))),
      saveInput_(cfg.getParameter<bool>("SaveInput")),
      model_(std::make_unique<L1TEmbeddingModelServer>(resolveModelName(cfg.getParameter<std::string>("ModelName")))) {
  produces<std::vector<float>>("L1TEmbeddingModel");
  if (saveInput_)
    produces<std::vector<float>>("L1TEmbeddingModelInput");
}

void L1TEmbeddingModelProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("PuppiCandidates", edm::InputTag("l1tLayer1Extended", "Puppi"));
  desc.add<std::string>("ModelName", "L1TEmbeddingModel");
  desc.add<bool>("SaveInput", false);
  descriptions.addWithDefaultLabel(desc);
}

void L1TEmbeddingModelProducer::produce(edm::StreamID, edm::Event &iEvent, const edm::EventSetup &) const {
  edm::Handle<std::vector<l1t::PFCandidate>> puppi;
  iEvent.getByToken(puppiToken_, puppi);

  std::vector<l1ct::PuppiObjEmu> particles;
  particles.reserve(puppi->size());
  for (const auto &cand : *puppi) {
    l1ct::PuppiObjEmu p;
    p.initFromBits(cand.encodedPuppi64());
    p.srcCand = &cand;
    particles.push_back(p);
  }

  std::vector<L1TEmbeddingModelServer::input_t> input;
  if (saveInput_)
    input.resize(L1TEmbeddingModelServer::N_inp);

  l1ct::L1TEmbeddingModelOutput output;
  l1ct::L1TEmbeddingModel_emu(particles, output, *model_, saveInput_ ? input.data() : nullptr);

  if (iEvent.id().event() == 1)
    std::cout << "L1TEmbeddingModelProducer event=" << iEvent.id().event() << " nPuppi=" << particles.size()
              << " embedding dimension=" << l1ct::L1TEmbeddingModelOutput::N_OUT << std::endl;

  std::vector<float> result;
  result.reserve(l1ct::L1TEmbeddingModelOutput::N_OUT);
  for (size_t i = 0; i < l1ct::L1TEmbeddingModelOutput::N_OUT; ++i)
    result.push_back(output.embedding[i].to_float());
  iEvent.put(std::make_unique<std::vector<float>>(std::move(result)), "L1TEmbeddingModel");

  if (saveInput_) {
    std::vector<float> inputOut;
    inputOut.reserve(L1TEmbeddingModelServer::N_inp);
    for (size_t i = 0; i < L1TEmbeddingModelServer::N_inp; ++i)
      inputOut.push_back(input[i].to_float());
    iEvent.put(std::make_unique<std::vector<float>>(std::move(inputOut)), "L1TEmbeddingModelInput");
  }
}

DEFINE_FWK_MODULE(L1TEmbeddingModelProducer);
