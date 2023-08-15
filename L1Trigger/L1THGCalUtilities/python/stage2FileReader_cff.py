import FWCore.ParameterSet.Config as cms

Stage2FileReader = cms.EDProducer('Stage2FileReader',
  files = cms.vstring("output_hgcAlgo.txt"),
  format = cms.untracked.string("EMPv2"),
  refClustersTag = cms.untracked.InputTag("l1tHGCalBackEndLayer2Producer", "HGCalBackendLayer2Processor3DClusteringSA"),
  sector = cms.untracked.uint32(0)
)
