import FWCore.ParameterSet.Config as cms

L1SeedConePFJetProducer = cms.EDProducer("L1SeedConePFJetProducer",
                           L1PFObjects = cms.InputTag("L1PFProducer","l1pfCandidates"),
                           nJets       = cms.uint32(10),
                           coneSize    = cms.double(0.4),
                           HW          = cms.bool(False),
                           debug       = cms.bool(False)
                         )

L1SeedConePFJetEmulatorProducer = cms.EDProducer("L1SeedConePFJetProducer",
                                   L1PFObjects = cms.InputTag("L1PFProducer","l1pfCandidates"),
                                   nJets       = cms.uint32(10),
                                   coneSize    = cms.double(0.4),
                                   HW          = cms.bool(True),
                                   debug       = cms.bool(False) 
                                  )
