import FWCore.ParameterSet.Config as cms

Stage2FileWriter = cms.EDAnalyzer('Stage2FileWriter',
  clusters = cms.untracked.InputTag("l1tHGCalBackEndLayer2Producer", "HGCalBackendLayer2Processor3DClusteringSA"),
  tmIndex = cms.untracked.uint32(0),
  format = cms.untracked.string("EMPv2")
)
