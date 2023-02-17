import FWCore.ParameterSet.Config as cms

Stage2FileWriter = cms.EDAnalyzer('Stage2FileWriter',
  clusters = cms.untracked.InputTag("l1tHGCalBackEndLayer2Producer", "HGCalBackendLayer2Processor3DClusteringSA"),
  outputFilename_clustersToL1T = cms.untracked.string("HGCS2OutputToL1TFile"),
  outputFilename_clusterSums = cms.untracked.string("HGCS2ClusterSumsInput"),
  format = cms.untracked.string("EMPv2")
)
