#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include <memory>
#include <string>
#include <vector>

class EmbeddingModelTableProducer : public edm::global::EDProducer<> {
public:
  explicit EmbeddingModelTableProducer(const edm::ParameterSet &);
  ~EmbeddingModelTableProducer() override = default;
  static void fillDescriptions(edm::ConfigurationDescriptions &);

private:
  void produce(edm::StreamID, edm::Event &, const edm::EventSetup &) const override;

  edm::EDGetTokenT<std::vector<float>> token_;
  std::string name_;
  std::string doc_;
};

EmbeddingModelTableProducer::EmbeddingModelTableProducer(const edm::ParameterSet &cfg)
    : token_(consumes<std::vector<float>>(cfg.getParameter<edm::InputTag>("src"))),
      name_(cfg.getParameter<std::string>("name")),
      doc_(cfg.getParameter<std::string>("doc")) {
  produces<nanoaod::FlatTable>();
}

void EmbeddingModelTableProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("l1tEmbeddingModelProducer", "L1TEmbeddingModel"));
  desc.add<std::string>("name", "L1Embedding");
  desc.add<std::string>("doc", "Embedding model output, 20-dim float vector");
  descriptions.addWithDefaultLabel(desc);
}

void EmbeddingModelTableProducer::produce(edm::StreamID, edm::Event &iEvent, const edm::EventSetup &) const {
  edm::Handle<std::vector<float>> handle;
  iEvent.getByToken(token_, handle);

  const auto &vec = *handle;
  const unsigned int n = vec.size();

  auto table = std::make_unique<nanoaod::FlatTable>(n, name_, false);
  table->setDoc(doc_);

  table->addColumn<float>("emb", vec, "Embedding output, 20-dim vector");

  iEvent.put(std::move(table));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(EmbeddingModelTableProducer);
