import FWCore.ParameterSet.Config as cms

# Specifications for Correlator Layer 1 output mapping onto links
barrelConfig_ = cms.PSet(
    partition = cms.string("Barrel"),
    nLinksPuppi = cms.uint32(6),
    nPuppiPerRegion = cms.uint32(18),
    nOutputFramesPerBX = cms.uint32(9),
    tmuxFactor = cms.uint32(6)
)

barrelPhiConfigs = [
    barrelConfig_.clone(
        outputRegions = cms.vuint32(*[3*ip+9*ie+i for ie in range(6) for i in range(3) ]),
        outputBoard = cms.int32(ip),
    ) for ip in range(3)
]

barrelPhiConfigs[2].outputBoard = 4

hgcalConfig_ = cms.PSet(
    partition = cms.string("HGCal"),
    nLinksPuppi = cms.uint32(3),
    nPuppiPerRegion = cms.uint32(18),
    nOutputFramesPerBX = cms.uint32(9),
    outputRegions = cms.vuint32(*[54 + i+9 for i in range(9)]),
    outputBoard = cms.int32(3),
    tmuxFactor = cms.uint32(6)
)

hgcalPosConfig = hgcalConfig_.clone(
    outputBoard = 3
)
hgcalNegConfig = hgcalConfig_.clone(
    outputRegions = [54 + i for i in range(9)],
    outputBoard = 2
)

hgcalNoTKConfig = cms.PSet(
    partition = cms.string("HGCalNoTk"),
    nLinksPuppi = cms.uint32(4),
    nPuppiPerRegion = cms.uint32(12),
    nOutputFramesPerBX = cms.uint32(9),
    outputRegions = cms.vuint32(*range(72,72+18)),
    outputBoard = cms.int32(5),
    tmuxFactor = cms.uint32(6)
)

hfConfig_ = cms.PSet(
    partition = cms.string("HF"),
    nLinksPuppi = cms.uint32(3),
    nPuppiPerRegion = cms.uint32(18),
    nOutputFramesPerBX = cms.uint32(9),
    tmuxFactor = cms.uint32(6)
)
hfConfigs = [
    hfConfig_.clone(
        outputRegions = cms.vuint32(*[90+9*ie+i for i in range(9)]),
        outputBoard = cms.int32(7 - ie),
    ) for ie in range(2)
]

linkConfigs = cms.VPSet(*barrelPhiConfigs, hgcalPosConfig, hgcalNegConfig, hgcalNoTKConfig, *hfConfigs)

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

################################################
# The following is an alternative configuration for the TM18 firmware, which has a different mapping of regions to output links.
################################################

barrelConfigTM18 = cms.PSet(
    partition = cms.string("Barrel"),
    nLinksPuppi = cms.uint32(6),
    nPuppiPerRegion = cms.uint32(18),
    nOutputFramesPerBX = cms.uint32(9),
    outputBoard = cms.int32(0),
    outputRegions = cms.vuint32(*range(54)),
    tmuxFactor = cms.uint32(18)
)

hgcalConfigTM18 = cms.PSet(
    partition = cms.string("HGCal"),
    nLinksPuppi = cms.uint32(2),
    nPuppiPerRegion = cms.uint32(18),
    nOutputFramesPerBX = cms.uint32(9),
    outputRegions = cms.vuint32(*range(54,72)),
    outputBoard = cms.int32(1),
    tmuxFactor = cms.uint32(18)
)

hgcalNoTKConfigTM18 = cms.PSet(
    partition = cms.string("HGCalNoTk"),
    nLinksPuppi = cms.uint32(2),
    nPuppiPerRegion = cms.uint32(12),
    nOutputFramesPerBX = cms.uint32(9),
    outputRegions = cms.vuint32(*range(72,72+18)),
    outputBoard = cms.int32(2),
    tmuxFactor = cms.uint32(18)
)

# hfConfigTM18 = ... no HF TM18 config for now
hfConfigTM18 = cms.PSet(
    partition = cms.string("HF"),
    nLinksPuppi = cms.uint32(1),
    nPuppiPerRegion = cms.uint32(6),
    nOutputFramesPerBX = cms.uint32(9),
    tmuxFactor = cms.uint32(6)
)

hfConfigTM18s = [
    hfConfigTM18.clone(
        outputRegions = cms.vuint32(*[90+9*ie+i for i in range(9)]),
        outputBoard = cms.int32(3 + ie),
    ) for ie in range(2)
]

linkConfigsTM18 = cms.VPSet(barrelConfigTM18, hgcalConfigTM18, hgcalNoTKConfigTM18, *hfConfigTM18s)
l1tDeregionizerProducerTM18 = l1tDeregionizerProducer.clone(linkConfigs = linkConfigsTM18)