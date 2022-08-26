import FWCore.ParameterSet.Config as cms

from L1Trigger.Phase2L1ParticleFlow.l1pfJetMet_cff import l1PFJetsExtendedTask

l1BJetProducerPuppi = cms.EDProducer("L1BJetProducer",
    jets = cms.InputTag("scPFL1PuppiExtended", ""),
    useRawPt = cms.bool(True),
    NNFileName = cms.string("L1Trigger/Phase2L1ParticleFlow/data/modelTT_PUP_Off_dXY_XYCut_Graph.pb"),
    maxJets = cms.int32(6),
    nParticles = cms.int32(10),
    minPt = cms.double(10),
    maxEta = cms.double(2.4),
    vtx = cms.InputTag("L1VertexFinderEmulator","l1verticesEmulation"),
    vtxEmulation = cms.bool(True),
    L1PFObjects = cms.InputTag("l1ctLayer1Extended","Puppi"),
)
l1BJetProducerPuppiCorrectedEmulator = l1BJetProducerPuppi.clone(
    jets = cms.InputTag("scPFL1PuppiExtendedCorrectedEmulator", ""),
    L1PFObjects = cms.InputTag("l1ctLayer2DeregionizerExtended","Puppi"),
)

l1BJetsTask = cms.Task(
    l1PFJetsExtendedTask, l1BJetProducerPuppi, l1BJetProducerPuppiCorrectedEmulator
)
