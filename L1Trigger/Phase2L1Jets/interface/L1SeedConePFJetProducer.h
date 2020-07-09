#ifndef L1Trigger_Phase2L1Jets_L1SeedConePFJetProducer_h
#define L1Trigger_Phase2L1Jets_L1SeedConePFJetProducer_h

#include <vector>
#include <numeric>

////////////////////
// FRAMEWORK HEADERS
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "DataFormats/L1TParticleFlow/interface/PFJet.h"
#include "DataFormats/Math/interface/deltaR.h"

class L1SeedConePFJetProducer : public edm::global::EDProducer<> {
   public:
       explicit L1SeedConePFJetProducer(const edm::ParameterSet&);
       ~L1SeedConePFJetProducer() override;

   private:

  /// ///////////////// ///
  /// MANDATORY METHODS ///
  void produce( edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup ) const override;
  /// ///////////////// ///

  float _coneSize;
  unsigned _nJets;
  edm::EDGetTokenT<std::vector<l1t::PFCandidate>> _l1PFToken;

  static l1t::PFJet makeJet(const std::vector< edm::Ptr< l1t::PFCandidate >> & parts) ;

};


#endif

