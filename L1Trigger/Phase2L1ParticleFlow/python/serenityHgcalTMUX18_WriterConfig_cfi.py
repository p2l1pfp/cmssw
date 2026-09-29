import FWCore.ParameterSet.Config as cms

writerConfig = cms.PSet(
    eventsPerFile = cms.uint32(4),
    fileFormat = cms.string('EMPv2'),
gmtLink = cms.int32(90),  # TM1:GMT -> 4*22+2
gmtNumberOfMuons = cms.uint32(12),
    gttLatency = cms.uint32(167),
gttLink = cms.int32(91),  # TM1:GTT -> 4*22+3
gttNumberOfPVs = cms.uint32(10),
    hgcSectors = cms.VPSet(
cms.PSet(
            hgcLinks = cms.vint32(75, 72, 68, 73)  # HGC_W1_N_L1:4*18+3, HGC_W1_N_L2:4*18+0, HGC_W1_N_L3:4*17+0, HGC_W1_N_L4:4*18+1
        ),
cms.PSet(
            hgcLinks = cms.vint32(69, 74, 70, 65)  # HGC_W2_N_L1:4*17+1, HGC_W2_N_L2:4*18+2, HGC_W2_N_L3:4*17+2, HGC_W2_N_L4:4*16+1
        ),
cms.PSet(
            hgcLinks = cms.vint32(71, 66, 64, 67)  # HGC_W3_N_L1:4*17+3, HGC_W3_N_L2:4*16+2, HGC_W3_N_L3:4*16+0, HGC_W3_N_L4:4*16+3
        ),
cms.PSet(
            hgcLinks = cms.vint32(60, 63, 59, 62)  # HGC_W1_P_L1:4*15+0, HGC_W1_P_L2:4*15+3, HGC_W1_P_L3:4*14+3, HGC_W1_P_L4:4*15+2
        ),
cms.PSet(
            hgcLinks = cms.vint32(58, 61, 57, 54)  # HGC_W2_P_L1:4*14+2, HGC_W2_P_L2:4*15+1, HGC_W2_P_L3:4*14+1, HGC_W2_P_L4:4*13+2
        ),
cms.PSet(
            hgcLinks = cms.vint32(56, 53, 55, 52)  # HGC_W3_P_L1:4*14+0, HGC_W3_P_L2:4*13+1, HGC_W3_P_L3:4*13+3, HGC_W3_P_L4:4*13+0
        )
),
    inputFileExtension = cms.string('txt.gz'),
    inputFileName = cms.string('l1HGCalTM18-inputs-vu13p'),
    maxLinesPerInputFile = cms.uint32(1191),
    maxLinesPerOutputFile = cms.uint32(1024),
    nEgammaObjectsOut = cms.uint32(16),
    nInputFramesPerBX = cms.uint32(9),
    nOutputFramesPerBX = cms.uint32(9),
    outputFileExtension = cms.string('txt.gz'),
    outputLinkEgamma = cms.int32(3),
    outputLinksPuppi = cms.vuint32(0, 1, 2),
    partition = cms.string('HGCal'),
    tfSectors = cms.VPSet(
cms.PSet(
            tfLink = cms.int32(7)  # TM1:TK_W1_N -> 4*1+3
        ),
cms.PSet(
            tfLink = cms.int32(4)  # TM1:TK_W2_N -> 4*1+0
        ),
cms.PSet(
            tfLink = cms.int32(8)  # TM1:TK_W3_N -> 4*2+0
        ),
cms.PSet(
            tfLink = cms.int32(5)  # TM1:TK_W4_N -> 4*1+1
        ),
cms.PSet(
            tfLink = cms.int32(9)  # TM1:TK_W5_N -> 4*2+1
        ),
cms.PSet(
            tfLink = cms.int32(6)  # TM1:TK_W6_N -> 4*1+2
        ),
cms.PSet(
            tfLink = cms.int32(10)  # TM1:TK_W7_N -> 4*2+2
        ),
cms.PSet(
            tfLink = cms.int32(13)  # TM1:TK_W8_N -> 4*3+1
        ),
cms.PSet(
            tfLink = cms.int32(11)  # TM1:TK_W9_N -> 4*2+3
        ),
cms.PSet(
            tfLink = cms.int32(123)  # TM1:TK_W1_P -> 4*30+3
        ),
cms.PSet(
            tfLink = cms.int32(120)  # TM1:TK_W2_P -> 4*30+0
        ),
cms.PSet(
            tfLink = cms.int32(116)  # TM1:TK_W3_P -> 4*29+0
        ),
cms.PSet(
            tfLink = cms.int32(121)  # TM1:TK_W4_P -> 4*30+1
        ),
cms.PSet(
            tfLink = cms.int32(117)  # TM1:TK_W5_P -> 4*29+1
        ),
cms.PSet(
            tfLink = cms.int32(122)  # TM1:TK_W6_P -> 4*30+2
        ),
cms.PSet(
            tfLink = cms.int32(118)  # TM1:TK_W7_P -> 4*29+2
        ),
cms.PSet(
            tfLink = cms.int32(113)  # TM1:TK_W8_P -> 4*28+1
        ),
cms.PSet(
            tfLink = cms.int32(119)  # TM1:TK_W9_P -> 4*29+3
        )
),
    tmuxFactor = cms.uint32(18)
)
