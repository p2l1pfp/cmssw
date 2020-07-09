#ifndef L1SCPFJET_PRDC_H
#define L1SCPFJET_PRDC_H

#include <vector>
#include <numeric>

////////////////////
// FRAMEWORK HEADERS
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Phase2L1ParticleFlow/interface/PFCandidate.h"
#include "DataFormats/Phase2L1ParticleFlow/interface/PFJet.h"
#include "DataFormats/Math/interface/deltaR.h"

class L1SeedConePFJetProducer : public edm::EDProducer {
   public:
  explicit L1SeedConePFJetProducer(const edm::ParameterSet&);

  ~L1SeedConePFJetProducer();

   private:

  /// ///////////////// ///
  /// MANDATORY METHODS ///
  virtual void beginRun( const edm::Run& run, const edm::EventSetup& iSetup );
  virtual void endRun( const edm::Run& run, const edm::EventSetup& iSetup );
  virtual void produce( edm::Event& iEvent, const edm::EventSetup& iSetup ) override;
  /// ///////////////// ///

  float _coneSize;
  unsigned _nJets;
  edm::EDGetTokenT<vector<l1t::PFCandidate>>_l1PFToken;

};


#endif

