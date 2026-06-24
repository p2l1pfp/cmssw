import FWCore.ParameterSet.Config as cms

# Specifications for Correlator Layer 1 output mapping onto links
barrelConfig = cms.PSet(
    partition = cms.string("Barrel"),
    nLinksPuppi = cms.uint32(6),
    nPuppiPerRegion = cms.uint32(18),
    nOutputFramesPerBX = cms.uint32(9),
    outputBoard = cms.int32(0),
    outputRegions = cms.vuint32(*range(54)),
    tmuxFactor = cms.uint32(18)
)

hgcalConfig = cms.PSet(
    partition = cms.string("HGCal"),
    nLinksPuppi = cms.uint32(2),
    nPuppiPerRegion = cms.uint32(18),
    nOutputFramesPerBX = cms.uint32(9),
    outputRegions = cms.vuint32(*range(54,72)),
    outputBoard = cms.int32(1),
    tmuxFactor = cms.uint32(18)
)

hgcalNoTKConfig = cms.PSet(
    partition = cms.string("HGCalNoTk"),
    nLinksPuppi = cms.uint32(2),
    nPuppiPerRegion = cms.uint32(12),
    nOutputFramesPerBX = cms.uint32(9),
    outputRegions = cms.vuint32(*range(72,72+18)),
    outputBoard = cms.int32(2),
    tmuxFactor = cms.uint32(18)
)

# HF configuration is best estimate at the moment
hfConfig = cms.PSet(
    partition = cms.string("HF"),
    nLinksPuppi = cms.uint32(1),
    nPuppiPerRegion = cms.uint32(6),
    nOutputFramesPerBX = cms.uint32(9),
    tmuxFactor = cms.uint32(6)
)

hfConfigs = [
    hfConfig.clone(
        outputRegions = cms.vuint32(*[90+9*ie+i for i in range(9)]),
        outputBoard = cms.int32(3 + ie),
    ) for ie in range(2)
]

linkConfigs = cms.VPSet(barrelConfig, hgcalConfig, hgcalNoTKConfig, *hfConfigs)


l1tDeregionizerProducer = cms.EDProducer("DeregionizerProducer",
                           RegionalPuppiCands  = cms.InputTag("l1tLayer1","PuppiRegional"),
                           nPuppiFinalBuffer   = cms.uint32(128),
                           nPuppiPerClk        = cms.uint32(6),
                           nPuppiFirstBuffers  = cms.uint32(12),
                           nPuppiSecondBuffers = cms.uint32(32),
                           nPuppiThirdBuffers  = cms.uint32(64),
                           nInputFramesPerBX   = cms.uint32(9),
                           linkConfigs         = linkConfigs,
                           writeInputPatternFiles = cms.bool(False),
                           inputPatternFilePSet = cms.PSet(
                                            gapLengthOutput = cms.uint32(0),
                                            TMUX = cms.uint32(6),
                                            maxLinesPerFile = cms.uint32(1024),
                                            outputFilename = cms.string("L1DeregionizerInput"),
                                            format = cms.string("EMPv2"),
                                            outputFileExtension = cms.string("txt.gz")
                                                     ),
                                                     inputPatternTimeSlices = cms.VPSet()
                         )

l1tDeregionizerProducerExtended = l1tDeregionizerProducer.clone(RegionalPuppiCands  = cms.InputTag("l1tLayer1Extended","PuppiRegional"))