import FWCore.ParameterSet.Config as cms

import math

from L1Trigger.Phase2L1ParticleFlow.l1tPFTracksFromL1Tracks_cfi import l1tPFTracksFromL1Tracks, l1tPFTracksFromL1TracksExtended
from L1Trigger.Phase2L1ParticleFlow.l1tPFClustersFromL1EGClusters_cfi import l1tPFClustersFromL1EGClusters
from L1Trigger.Phase2L1ParticleFlow.pfClustersFromCombinedCalo_cff import l1tPFClustersFromCombinedCaloHCal, l1tPFClustersFromCombinedCaloHF
from L1Trigger.Phase2L1ParticleFlow.l1tPFClustersFromHGC3DClusters_cfi import l1tPFClustersFromHGC3DClusters

from L1Trigger.Phase2L1ParticleFlow.l1TkEgAlgoEmulator_cfi import tkEgAlgoParameters,tkEgSorterParameters
from L1Trigger.VertexFinder.l1tVertexProducer_cfi import l1tVertexProducer

l1tLayer1Barrel = cms.EDProducer("L1TCorrelatorLayer1Producer",
    tracks = cms.InputTag('l1tPFTracksFromL1Tracks'),
    muons = cms.InputTag('l1tSAMuonsGmt','prompt'),
    emClusters = cms.VInputTag(cms.InputTag('l1tPFClustersFromL1EGClusters:selected')),
    hadClusters = cms.VInputTag(cms.InputTag('l1tPFClustersFromCombinedCaloHCal:calibrated')),
    vtxCollection = cms.InputTag("l1tVertexFinderEmulator","L1VerticesEmulation"),
    nVtx = cms.int32(1),    
    emPtCut = cms.double(0.5),
    hadPtCut = cms.double(1.0),
    trkPtCut = cms.double(2.0),
    trackInputConversionAlgo = cms.string("Emulator"),
    trackInputConversionParameters = cms.PSet(
        region = cms.string("barrel"),
        slimDataFormat = cms.bool(True),
        ptLUTBits = cms.uint32(11),
        etaLUTBits = cms.uint32(10),
        etaShift = cms.uint32(15-10),
        phiBits = cms.uint32(10),
        z0Bits = cms.uint32(12),
        dEtaBarrelBits = cms.uint32(8),
        dEtaBarrelZ0PreShift = cms.uint32(2),
        dEtaBarrelZ0PostShift = cms.uint32(2),
        dPhiBarrelBits = cms.uint32(4),
        dPhiBarrelRInvPreShift = cms.uint32(4),
        dPhiBarrelRInvPostShift = cms.uint32(4),
        ),
    muonInputConversionAlgo = cms.string("Emulator"),
    regionizerAlgo = cms.string("Ideal"),
    pfAlgo = cms.string("PFAlgo3"),
    pfAlgoParameters = cms.PSet(
        nTrack = cms.uint32(25), 
        nCalo = cms.uint32(18), 
        nMu = cms.uint32(2), 
        nSelCalo = cms.uint32(18), 
        nEmCalo = cms.uint32(12), 
        nPhoton = cms.uint32(12), 
        nAllNeutral = cms.uint32(25), 
        trackMuDR    = cms.double(0.2), # accounts for poor resolution of standalone, and missing propagations
        trackEmDR   = cms.double(0.04), # 1 Ecal crystal size is 0.02, and ~2 cm in HGCal is ~0.007
        emCaloDR    = cms.double(0.10),    # 1 Hcal tower size is ~0.09
        trackCaloDR = cms.double(0.15),
        maxInvisiblePt = cms.double(10.0), # max allowed pt of a track with no calo energy
        tightTrackMaxInvisiblePt = cms.double(20),
        caloResolution = cms.PSet(
            etaBins = cms.vdouble( 0.700,  1.200,  1.600),
            offset  = cms.vdouble( 2.909,  2.864,  0.294),
            scale   = cms.vdouble( 0.119,  0.127,  0.442),
        ),
    ),
    puAlgo = cms.string("LinearizedPuppi"),
    puAlgoParameters = cms.PSet(
        nTrack = cms.uint32(22), 
        nIn = cms.uint32(25), 
        nOut = cms.uint32(25), 
        nVtx = cms.uint32(1),
        nFinalSort = cms.uint32(18), 
        finalSortAlgo = cms.string("Insertion"),
        dZ     = cms.double(0.5),
        dr     = cms.double(0.3),
        drMin  = cms.double(0.07),
        ptMax  = cms.double(50.),
        absEtaCuts         = cms.vdouble( ), # just one bin, so no edge needd
        ptCut             = cms.vdouble( 1.0 ),
        ptSlopes           = cms.vdouble( 0.3 ), # coefficient for pT
        ptSlopesPhoton    = cms.vdouble( 0.3 ),
        ptZeros            = cms.vdouble( 4.0 ), # ballpark pT from PU
        ptZerosPhoton     = cms.vdouble( 2.5 ),
        alphaSlopes        = cms.vdouble( 0.7 ), # coefficient for alpha
        alphaZeros         = cms.vdouble( 6.0 ), # ballpark alpha from PU
        alphaCrop         = cms.vdouble(  4  ), # max. absolute value for alpha term
        priors             = cms.vdouble( 5.0 ),
        priorsPhoton      = cms.vdouble( 1.0 ),
        useAssociationNetwork = cms.bool(True), #Enable Association Network
        associationThreshold = cms.double(0.1), #Association Network threshold for PV tracks
        associationGraph = cms.FileInPath("L1Trigger/L1TTrackMatch/data/NNVtx_AssociationModelGraph.pb"),
        associationNetworkZ0binning = l1tVertexProducer.VertexReconstruction.FH_HistogramParameters, #Z0 binning used for setting the input feature digitisation
        associationNetworkEtaBounds = cms.vdouble(0.0, 0.01984126984126984, 0.03968253968253968, 0.05952380952380952, 0.07936507936507936, 0.0992063492063492, 0.11904761904761904, 0.1388888888888889, 0.15873015873015872, 0.17857142857142855, 0.1984126984126984, 0.21825396825396826, 0.23809523809523808, 0.2579365079365079, 0.2777777777777778, 0.2976190476190476, 0.31746031746031744, 0.33730158730158727, 0.3571428571428571, 0.376984126984127, 0.3968253968253968, 0.41666666666666663, 0.4365079365079365, 0.45634920634920634, 0.47619047619047616, 0.496031746031746, 0.5158730158730158, 0.5357142857142857, 0.5555555555555556, 0.5753968253968254, 0.5952380952380952, 0.615079365079365, 0.6349206349206349, 0.6547619047619048, 0.6746031746031745, 0.6944444444444444, 0.7142857142857142, 0.7341269841269841, 0.753968253968254, 0.7738095238095237, 0.7936507936507936, 0.8134920634920635, 0.8333333333333333, 0.8531746031746031, 0.873015873015873, 0.8928571428571428, 0.9126984126984127, 0.9325396825396824, 0.9523809523809523, 0.9722222222222222, 0.992063492063492, 1.0119047619047619, 1.0317460317460316, 1.0515873015873016, 1.0714285714285714, 1.0912698412698412, 1.1111111111111112, 1.130952380952381, 1.1507936507936507, 1.1706349206349205, 1.1904761904761905, 1.2103174603174602, 1.23015873015873, 1.25, 1.2698412698412698, 1.2896825396825395, 1.3095238095238095, 1.3293650793650793, 1.349206349206349, 1.369047619047619, 1.3888888888888888, 1.4087301587301586, 1.4285714285714284, 1.4484126984126984, 1.4682539682539681, 1.488095238095238, 1.507936507936508, 1.5277777777777777, 1.5476190476190474, 1.5674603174603174, 1.5873015873015872, 1.607142857142857, 1.626984126984127, 1.6468253968253967, 1.6666666666666665, 1.6865079365079365, 1.7063492063492063, 1.726190476190476, 1.746031746031746, 1.7658730158730158, 1.7857142857142856, 1.8055555555555554, 1.8253968253968254, 1.8452380952380951, 1.865079365079365, 1.8849206349206349, 1.9047619047619047, 1.9246031746031744, 1.9444444444444444, 1.9642857142857142, 1.984126984126984, 2.003968253968254, 2.0238095238095237, 2.0436507936507935, 2.0634920634920633, 2.083333333333333, 2.1031746031746033, 2.123015873015873, 2.142857142857143, 2.1626984126984126, 2.1825396825396823, 2.202380952380952, 2.2222222222222223, 2.242063492063492, 2.261904761904762, 2.2817460317460316, 2.3015873015873014, 2.321428571428571, 2.341269841269841, 2.361111111111111, 2.380952380952381, 2.4007936507936507, 2.4206349206349205, 2.4404761904761902, 2.46031746031746, 2.4801587301587302, 2.5), #Eta bounds used to set z0 resolution input feature
        associationNetworkZ0ResBins = cms.vdouble(127.0, 126.0, 126.0, 126.0, 125.0, 124.0, 123.0, 122.0, 120.0, 119.0, 117.0, 115.0, 114.0, 112.0, 110.0, 107.0, 105.0, 103.0, 101.0, 98.0, 96.0, 94.0, 91.0, 89.0, 87.0, 85.0, 82.0, 80.0, 78.0, 76.0, 74.0, 72.0, 70.0, 68.0, 66.0, 64.0, 62.0, 61.0, 59.0, 57.0, 56.0, 54.0, 53.0, 51.0, 50.0, 48.0, 47.0, 46.0, 45.0, 43.0, 42.0, 41.0, 40.0, 39.0, 38.0, 37.0, 36.0, 35.0, 34.0, 33.0, 33.0, 32.0, 31.0, 30.0, 30.0, 29.0, 28.0, 28.0, 27.0, 26.0, 26.0, 25.0, 24.0, 24.0, 23.0, 23.0, 22.0, 22.0, 21.0, 21.0, 21.0, 20.0, 20.0, 19.0, 19.0, 18.0, 18.0, 18.0, 17.0, 17.0, 17.0, 16.0, 16.0, 16.0, 15.0, 15.0, 15.0, 15.0, 14.0, 14.0, 14.0, 14.0, 13.0, 13.0, 13.0, 13.0, 12.0, 12.0, 12.0, 12.0, 12.0, 11.0, 11.0, 11.0, 11.0, 11.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 9.0, 9.0, 9.0, 9.0, 9.0, 0.0), #z0 resolution input feature bins
    ),
    tkEgAlgoParameters = tkEgAlgoParameters.clone(
        nTRACK = 25,
        nTRACK_EGIN = 13,
        nEMCALO_EGIN = 10,
        nEM_EGOUT = 10,
    ),
    tkEgSorterAlgo = cms.string("Barrel"),
    tkEgSorterParameters = tkEgSorterParameters.clone(
        nObjToSort = 10
    ),
    caloSectors = cms.VPSet(
        cms.PSet( 
            etaBoundaries = cms.vdouble(-1.5, 1.5),
            phiSlices     = cms.uint32(3),
        )
    ),
    regions = cms.VPSet(
        cms.PSet( 
            etaBoundaries = cms.vdouble(-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5),
            phiSlices     = cms.uint32(9),
        ),
    ),
    boards = cms.VPSet(
        cms.PSet(
              regions = cms.vuint32(*[0+9*ie+i for ie in range(6) for i in range(3)])), # phi splitting
        cms.PSet(
              regions = cms.vuint32(*[3+9*ie+i for ie in range(6) for i in range(3)])), # phi splitting
        cms.PSet(
              regions = cms.vuint32(*[6+9*ie+i for ie in range(6) for i in range(3)])), # phi splitting
    )
)

l1tLayer1BarrelExtended = l1tLayer1Barrel.clone(tracks = cms.InputTag('l1tPFTracksFromL1TracksExtended'))

_hgcalSectors = cms.VPSet(
    cms.PSet( 
        etaBoundaries = cms.vdouble(-3.0, -1.5),
        phiSlices     = cms.uint32(3),
        phiZero       = cms.double(math.pi/2) # The edge of the 0th HGCal sectors is at 30 deg, the center at 30+120/2 = 90 = pi/2
    ),
    cms.PSet( 
        etaBoundaries = cms.vdouble(+1.5, +3.0),
        phiSlices     = cms.uint32(3),
        phiZero       = cms.double(math.pi/2) # As above
    )

)

l1tLayer1HGCal = cms.EDProducer("L1TCorrelatorLayer1Producer",
    tracks = cms.InputTag('l1tPFTracksFromL1Tracks'),
    muons = cms.InputTag('l1tSAMuonsGmt','prompt'),
    emClusters = cms.VInputTag(cms.InputTag('l1tPFClustersFromHGC3DClusters:egamma')), # used only for E/gamma
    hadClusters = cms.VInputTag(cms.InputTag('l1tPFClustersFromHGC3DClusters')),
    vtxCollection = cms.InputTag("l1tVertexFinderEmulator","L1VerticesEmulation"),
    nVtx = cms.int32(1),    
    emPtCut = cms.double(0.5),
    hadPtCut = cms.double(1.0),
    trkPtCut = cms.double(2.0),
    trackInputConversionAlgo = cms.string("Emulator"),
    trackInputConversionParameters = cms.PSet(
        region = cms.string("endcap"),
        slimDataFormat = cms.bool(False),
        ptLUTBits = cms.uint32(11),
        etaLUTBits = cms.uint32(11),
        etaShift = cms.uint32(15-11),
        etaPostOffs = cms.int32(150),
        phiBits = cms.uint32(10),
        z0Bits = cms.uint32(12),
        dEtaHGCalBits = cms.uint32(10),
        dEtaHGCalZ0PreShift = cms.uint32(2),
        dEtaHGCalRInvPreShift = cms.uint32(6),
        dEtaHGCalLUTBits = cms.uint32(10),
        dEtaHGCalLUTShift = cms.uint32(2),
        dPhiHGCalBits = cms.uint32(4),
        dPhiHGCalZ0PreShift = cms.uint32(4),
        dPhiHGCalZ0PostShift = cms.uint32(6),
        dPhiHGCalRInvShift = cms.uint32(4),
        dPhiHGCalTanlInvShift = cms.uint32(22),
        dPhiHGCalTanlLUTBits = cms.uint32(10),
        ),
    muonInputConversionAlgo = cms.string("Emulator"),
    hgcalInputConversionAlgo = cms.string("Emulator"),
    regionizerAlgo = cms.string("Multififo"),
    regionizerAlgoParameters = cms.PSet(
        useAlsoVtxCoords = cms.bool(True),
        nEndcaps = cms.uint32(2),
        nClocks = cms.uint32(54),
        nTkLinks = cms.uint32(2),
        nCaloLinks = cms.uint32(3),
        nTrack = cms.uint32(30),
        nCalo = cms.uint32(20),
        nEmCalo = cms.uint32(10),
        nMu = cms.uint32(4),
        egInterceptMode = cms.PSet(
            afterFifo = cms.bool(True),
            emIDMask = cms.uint32(0x1E),
            nHADCALO_IN = cms.uint32(20),
            nEMCALO_OUT = cms.uint32(10),
            )
        ),
    pfAlgo = cms.string("PFAlgo2HGC"),
    pfAlgoParameters = cms.PSet(
        nTrack = cms.uint32(30),
        nCalo = cms.uint32(20),
        nMu = cms.uint32(4),
        nSelCalo = cms.uint32(20),
        trackMuDR    = cms.double(0.2), # accounts for poor resolution of standalone, and missing propagations
        trackCaloDR = cms.double(0.1),
        maxInvisiblePt = cms.double(10.0), # max allowed pt of a track with no calo energy
        tightTrackMaxInvisiblePt = cms.double(20),
        caloResolution = cms.PSet(
            etaBins = cms.vdouble( 1.700,  1.900,  2.200,  2.500,  2.800,  2.900),
            offset  = cms.vdouble( 1.793,  1.827,  2.363,  2.538,  2.812,  2.642),
            scale   = cms.vdouble( 0.138,  0.137,  0.124,  0.115,  0.106,  0.121),
        ),
    ),
    puAlgo = cms.string("LinearizedPuppi"),
    puAlgoParameters = cms.PSet(
        nTrack = cms.uint32(30),
        nIn = cms.uint32(20),
        nOut = cms.uint32(20),
        nVtx        = cms.uint32(1),    
        nFinalSort = cms.uint32(18), 
        finalSortAlgo = cms.string("FoldedHybrid"),
        dZ     = cms.double(1.33),
        dr     = cms.double(0.3),
        drMin  = cms.double(0.04),
        ptMax  = cms.double(50.),
        absEtaCuts         = cms.vdouble( 2.0 ), # two bins in the tracker (different eta); give only the one boundary between them 
        ptCut             = cms.vdouble( 1.0, 2.0 ),
        ptSlopes           = cms.vdouble( 0.3, 0.3 ), # coefficient for pT
        ptSlopesPhoton    = cms.vdouble( 0.4, 0.4 ), #When e/g ID not applied, use: cms.vdouble( 0.3, 0.3, 0.3 ),
        ptZeros            = cms.vdouble( 5.0, 7.0 ), # ballpark pT from PU
        ptZerosPhoton     = cms.vdouble( 3.0, 4.0 ),
        alphaSlopes        = cms.vdouble( 1.5, 1.5 ),
        alphaZeros         = cms.vdouble( 6.0, 6.0 ),
        alphaCrop         = cms.vdouble(  3 ,  3  ), # max. absolute value for alpha term
        priors             = cms.vdouble( 5.0, 5.0 ),
        priorsPhoton      = cms.vdouble( 1.5, 1.5 ), #When e/g ID not applied, use: cms.vdouble( 3.5, 3.5, 7.0 ),
        useAssociationNetwork = cms.bool(True), #Enable Association Network
        associationThreshold = cms.double(0.1), #Association Network threshold for PV tracks
        associationGraph = cms.FileInPath("L1Trigger/L1TTrackMatch/data/NNVtx_AssociationModelGraph.pb"),
        associationNetworkZ0binning = l1tVertexProducer.VertexReconstruction.FH_HistogramParameters, #Z0 binning used for setting the input feature digitisation
        associationNetworkEtaBounds = cms.vdouble(0.0, 0.01984126984126984, 0.03968253968253968, 0.05952380952380952, 0.07936507936507936, 0.0992063492063492, 0.11904761904761904, 0.1388888888888889, 0.15873015873015872, 0.17857142857142855, 0.1984126984126984, 0.21825396825396826, 0.23809523809523808, 0.2579365079365079, 0.2777777777777778, 0.2976190476190476, 0.31746031746031744, 0.33730158730158727, 0.3571428571428571, 0.376984126984127, 0.3968253968253968, 0.41666666666666663, 0.4365079365079365, 0.45634920634920634, 0.47619047619047616, 0.496031746031746, 0.5158730158730158, 0.5357142857142857, 0.5555555555555556, 0.5753968253968254, 0.5952380952380952, 0.615079365079365, 0.6349206349206349, 0.6547619047619048, 0.6746031746031745, 0.6944444444444444, 0.7142857142857142, 0.7341269841269841, 0.753968253968254, 0.7738095238095237, 0.7936507936507936, 0.8134920634920635, 0.8333333333333333, 0.8531746031746031, 0.873015873015873, 0.8928571428571428, 0.9126984126984127, 0.9325396825396824, 0.9523809523809523, 0.9722222222222222, 0.992063492063492, 1.0119047619047619, 1.0317460317460316, 1.0515873015873016, 1.0714285714285714, 1.0912698412698412, 1.1111111111111112, 1.130952380952381, 1.1507936507936507, 1.1706349206349205, 1.1904761904761905, 1.2103174603174602, 1.23015873015873, 1.25, 1.2698412698412698, 1.2896825396825395, 1.3095238095238095, 1.3293650793650793, 1.349206349206349, 1.369047619047619, 1.3888888888888888, 1.4087301587301586, 1.4285714285714284, 1.4484126984126984, 1.4682539682539681, 1.488095238095238, 1.507936507936508, 1.5277777777777777, 1.5476190476190474, 1.5674603174603174, 1.5873015873015872, 1.607142857142857, 1.626984126984127, 1.6468253968253967, 1.6666666666666665, 1.6865079365079365, 1.7063492063492063, 1.726190476190476, 1.746031746031746, 1.7658730158730158, 1.7857142857142856, 1.8055555555555554, 1.8253968253968254, 1.8452380952380951, 1.865079365079365, 1.8849206349206349, 1.9047619047619047, 1.9246031746031744, 1.9444444444444444, 1.9642857142857142, 1.984126984126984, 2.003968253968254, 2.0238095238095237, 2.0436507936507935, 2.0634920634920633, 2.083333333333333, 2.1031746031746033, 2.123015873015873, 2.142857142857143, 2.1626984126984126, 2.1825396825396823, 2.202380952380952, 2.2222222222222223, 2.242063492063492, 2.261904761904762, 2.2817460317460316, 2.3015873015873014, 2.321428571428571, 2.341269841269841, 2.361111111111111, 2.380952380952381, 2.4007936507936507, 2.4206349206349205, 2.4404761904761902, 2.46031746031746, 2.4801587301587302, 2.5), #Eta bounds used to set z0 resolution input feature
        associationNetworkZ0ResBins = cms.vdouble(127.0, 126.0, 126.0, 126.0, 125.0, 124.0, 123.0, 122.0, 120.0, 119.0, 117.0, 115.0, 114.0, 112.0, 110.0, 107.0, 105.0, 103.0, 101.0, 98.0, 96.0, 94.0, 91.0, 89.0, 87.0, 85.0, 82.0, 80.0, 78.0, 76.0, 74.0, 72.0, 70.0, 68.0, 66.0, 64.0, 62.0, 61.0, 59.0, 57.0, 56.0, 54.0, 53.0, 51.0, 50.0, 48.0, 47.0, 46.0, 45.0, 43.0, 42.0, 41.0, 40.0, 39.0, 38.0, 37.0, 36.0, 35.0, 34.0, 33.0, 33.0, 32.0, 31.0, 30.0, 30.0, 29.0, 28.0, 28.0, 27.0, 26.0, 26.0, 25.0, 24.0, 24.0, 23.0, 23.0, 22.0, 22.0, 21.0, 21.0, 21.0, 20.0, 20.0, 19.0, 19.0, 18.0, 18.0, 18.0, 17.0, 17.0, 17.0, 16.0, 16.0, 16.0, 15.0, 15.0, 15.0, 15.0, 14.0, 14.0, 14.0, 14.0, 13.0, 13.0, 13.0, 13.0, 12.0, 12.0, 12.0, 12.0, 12.0, 11.0, 11.0, 11.0, 11.0, 11.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 9.0, 9.0, 9.0, 9.0, 9.0, 0.0), #z0 resolution input feature bins
    ),
    tkEgAlgoParameters = tkEgAlgoParameters.clone(
        nTRACK = 30,
        nTRACK_EGIN = 10,
        nEMCALO_EGIN = 10, 
        nEM_EGOUT = 5,
        doBremRecovery = True,
        doEndcapHwQual = True,
        writeBeforeBremRecovery = False,
        writeEGSta = True,
        doCompositeTkEle = True,
        trkQualityPtMin = 0.), # This should be 10 GeV when doCompositeTkEle = False
    tkEgSorterAlgo = cms.string("Endcap"),
    tkEgSorterParameters = tkEgSorterParameters.clone(
        nObjToSort = 5
    ),
    caloSectors = _hgcalSectors,
    regions = cms.VPSet(
        cms.PSet( 
            etaBoundaries = cms.vdouble(-2.5, -1.5),
            phiSlices     = cms.uint32(9),
        ),
        cms.PSet( 
            etaBoundaries = cms.vdouble(+1.5, +2.5),
            phiSlices     = cms.uint32(9),
        )

    ),
    boards = cms.VPSet(
        cms.PSet(
            regions = cms.vuint32(range(0, 9))),
        cms.PSet(
            regions = cms.vuint32(range(9, 18))),
    ),
    writeRawHgcalCluster = cms.untracked.bool(True)
)


l1tLayer1HGCalExtended = l1tLayer1HGCal.clone(tracks = ('l1tPFTracksFromL1TracksExtended'))

l1tLayer1HGCalElliptic = l1tLayer1HGCal.clone(
    tkEgAlgoParameters = l1tLayer1HGCal.tkEgAlgoParameters.clone(
        doCompositeTkEle = False,
        trkQualityPtMin = 10.)
)

l1tLayer1HGCalNoTK = cms.EDProducer("L1TCorrelatorLayer1Producer",
    muons = cms.InputTag('l1tSAMuonsGmt','prompt'),
    emClusters = cms.VInputTag(cms.InputTag('l1tPFClustersFromHGC3DClusters:egamma')), # used only for E/gamma
    hadClusters = cms.VInputTag(cms.InputTag('l1tPFClustersFromHGC3DClusters')),
    vtxCollection = cms.InputTag("l1tVertexFinderEmulator","L1VerticesEmulation"),
    nVtx = cms.int32(1),        
    emPtCut = cms.double(0.5),
    hadPtCut = cms.double(1.0),
    trkPtCut = cms.double(2.0),
    muonInputConversionAlgo = cms.string("Emulator"),
    hgcalInputConversionAlgo = cms.string("Emulator"),
    hgcalInputConversionParameters = cms.PSet(
        slim = cms.bool(True)
    ),
    regionizerAlgo = cms.string("Multififo"),
    regionizerAlgoParameters = cms.PSet(
        useAlsoVtxCoords = cms.bool(True),
        nEndcaps = cms.uint32(2),
        nClocks = cms.uint32(54),
        nTkLinks = cms.uint32(0),
        nCaloLinks = cms.uint32(3),
        nTrack = cms.uint32(0),
        nCalo = cms.uint32(12),
        nEmCalo = cms.uint32(12),
        nMu = cms.uint32(4),
        egInterceptMode = cms.PSet(
            afterFifo = cms.bool(True),
            emIDMask = cms.uint32(0x1E),
            nHADCALO_IN = cms.uint32(12),
            nEMCALO_OUT = cms.uint32(12),
            )
        ),
    pfAlgo = cms.string("PFAlgoDummy"),
    pfAlgoParameters = cms.PSet(
        nCalo = cms.uint32(12), 
        nMu = cms.uint32(4), # unused
    ),
    puAlgo = cms.string("LinearizedPuppi"),
    puAlgoParameters = cms.PSet(
        nTrack = cms.uint32(0),  # unused
        nIn = cms.uint32(12), 
        nOut = cms.uint32(12), 
        nFinalSort = cms.uint32(12), # to be tuned
        finalSortAlgo = cms.string("Hybrid"), 
        nVtx = cms.uint32(1),    
        dZ     = cms.double(1.33),
        dr     = cms.double(0.3),
        drMin  = cms.double(0.04),
        ptMax  = cms.double(50.),
        absEtaCuts         = cms.vdouble( ), # just one bin
        ptCut             = cms.vdouble( 4.0 ),
        ptSlopes           = cms.vdouble( 0.3 ), # coefficient for pT
        ptSlopesPhoton    = cms.vdouble( 0.4 ), #When e/g ID not applied, use: cms.vdouble( 0.3, 0.3, 0.3 ),
        ptZeros            = cms.vdouble( 9.0 ), # ballpark pT from PU
        ptZerosPhoton     = cms.vdouble( 5.0 ),
        alphaSlopes        = cms.vdouble( 2.2 ),
        alphaZeros         = cms.vdouble( 9.0 ),
        alphaCrop         = cms.vdouble(  4  ), # max. absolute value for alpha term
        priors             = cms.vdouble( 7.0 ),
        priorsPhoton      = cms.vdouble( 5.0 ), #When e/g ID not applied, use: cms.vdouble( 3.5, 3.5, 7.0 ),
        useAssociationNetwork = cms.bool(True), #Enable Association Network
        associationThreshold = cms.double(0.1), #Association Network threshold for PV tracks
        associationGraph = cms.FileInPath("L1Trigger/L1TTrackMatch/data/NNVtx_AssociationModelGraph.pb"),
        associationNetworkZ0binning = l1tVertexProducer.VertexReconstruction.FH_HistogramParameters, #Z0 binning used for setting the input feature digitisation
        associationNetworkEtaBounds = cms.vdouble(0.0, 0.01984126984126984, 0.03968253968253968, 0.05952380952380952, 0.07936507936507936, 0.0992063492063492, 0.11904761904761904, 0.1388888888888889, 0.15873015873015872, 0.17857142857142855, 0.1984126984126984, 0.21825396825396826, 0.23809523809523808, 0.2579365079365079, 0.2777777777777778, 0.2976190476190476, 0.31746031746031744, 0.33730158730158727, 0.3571428571428571, 0.376984126984127, 0.3968253968253968, 0.41666666666666663, 0.4365079365079365, 0.45634920634920634, 0.47619047619047616, 0.496031746031746, 0.5158730158730158, 0.5357142857142857, 0.5555555555555556, 0.5753968253968254, 0.5952380952380952, 0.615079365079365, 0.6349206349206349, 0.6547619047619048, 0.6746031746031745, 0.6944444444444444, 0.7142857142857142, 0.7341269841269841, 0.753968253968254, 0.7738095238095237, 0.7936507936507936, 0.8134920634920635, 0.8333333333333333, 0.8531746031746031, 0.873015873015873, 0.8928571428571428, 0.9126984126984127, 0.9325396825396824, 0.9523809523809523, 0.9722222222222222, 0.992063492063492, 1.0119047619047619, 1.0317460317460316, 1.0515873015873016, 1.0714285714285714, 1.0912698412698412, 1.1111111111111112, 1.130952380952381, 1.1507936507936507, 1.1706349206349205, 1.1904761904761905, 1.2103174603174602, 1.23015873015873, 1.25, 1.2698412698412698, 1.2896825396825395, 1.3095238095238095, 1.3293650793650793, 1.349206349206349, 1.369047619047619, 1.3888888888888888, 1.4087301587301586, 1.4285714285714284, 1.4484126984126984, 1.4682539682539681, 1.488095238095238, 1.507936507936508, 1.5277777777777777, 1.5476190476190474, 1.5674603174603174, 1.5873015873015872, 1.607142857142857, 1.626984126984127, 1.6468253968253967, 1.6666666666666665, 1.6865079365079365, 1.7063492063492063, 1.726190476190476, 1.746031746031746, 1.7658730158730158, 1.7857142857142856, 1.8055555555555554, 1.8253968253968254, 1.8452380952380951, 1.865079365079365, 1.8849206349206349, 1.9047619047619047, 1.9246031746031744, 1.9444444444444444, 1.9642857142857142, 1.984126984126984, 2.003968253968254, 2.0238095238095237, 2.0436507936507935, 2.0634920634920633, 2.083333333333333, 2.1031746031746033, 2.123015873015873, 2.142857142857143, 2.1626984126984126, 2.1825396825396823, 2.202380952380952, 2.2222222222222223, 2.242063492063492, 2.261904761904762, 2.2817460317460316, 2.3015873015873014, 2.321428571428571, 2.341269841269841, 2.361111111111111, 2.380952380952381, 2.4007936507936507, 2.4206349206349205, 2.4404761904761902, 2.46031746031746, 2.4801587301587302, 2.5), #Eta bounds used to set z0 resolution input feature
        associationNetworkZ0ResBins = cms.vdouble(127.0, 126.0, 126.0, 126.0, 125.0, 124.0, 123.0, 122.0, 120.0, 119.0, 117.0, 115.0, 114.0, 112.0, 110.0, 107.0, 105.0, 103.0, 101.0, 98.0, 96.0, 94.0, 91.0, 89.0, 87.0, 85.0, 82.0, 80.0, 78.0, 76.0, 74.0, 72.0, 70.0, 68.0, 66.0, 64.0, 62.0, 61.0, 59.0, 57.0, 56.0, 54.0, 53.0, 51.0, 50.0, 48.0, 47.0, 46.0, 45.0, 43.0, 42.0, 41.0, 40.0, 39.0, 38.0, 37.0, 36.0, 35.0, 34.0, 33.0, 33.0, 32.0, 31.0, 30.0, 30.0, 29.0, 28.0, 28.0, 27.0, 26.0, 26.0, 25.0, 24.0, 24.0, 23.0, 23.0, 22.0, 22.0, 21.0, 21.0, 21.0, 20.0, 20.0, 19.0, 19.0, 18.0, 18.0, 18.0, 17.0, 17.0, 17.0, 16.0, 16.0, 16.0, 15.0, 15.0, 15.0, 15.0, 14.0, 14.0, 14.0, 14.0, 13.0, 13.0, 13.0, 13.0, 12.0, 12.0, 12.0, 12.0, 12.0, 11.0, 11.0, 11.0, 11.0, 11.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 9.0, 9.0, 9.0, 9.0, 9.0, 0.0), #z0 resolution input feature bins
    ),
    tkEgAlgoParameters = tkEgAlgoParameters.clone(
        nTRACK = 30,
        nTRACK_EGIN = 10,
        nEMCALO_EGIN = 10, 
        nEM_EGOUT = 5,
        doBremRecovery = True,
        doEndcapHwQual = True,
        writeBeforeBremRecovery = False,
        writeEGSta = True),
    tkEgSorterAlgo = cms.string("Endcap"),
    tkEgSorterParameters = tkEgSorterParameters.clone(
        nObjToSort = 5
    ),
    caloSectors = _hgcalSectors,
    regions = cms.VPSet(
        cms.PSet( 
            etaBoundaries = cms.vdouble(-3.0, -2.5),
            phiSlices     = cms.uint32(9),
        ),
        cms.PSet( 
            etaBoundaries = cms.vdouble(+2.5, +3.0),
            phiSlices     = cms.uint32(9),
        )

    ),
    boards = cms.VPSet(
        cms.PSet(regions = cms.vuint32(range(0,18))),
    ),
    writeRawHgcalCluster = cms.untracked.bool(True)
)

l1tLayer1HF = cms.EDProducer("L1TCorrelatorLayer1Producer",
    muons = cms.InputTag('l1tSAMuonsGmt','prompt'),
    hadClusters = cms.VInputTag(cms.InputTag('l1tPFClustersFromCombinedCaloHF:calibrated')),
    vtxCollection = cms.InputTag("l1tVertexFinderEmulator","L1VerticesEmulation"),
    nVtx = cms.int32(1),    
    emPtCut = cms.double(0.5),
    hadPtCut = cms.double(15.0),
    trkPtCut = cms.double(2.0),
    pfAlgo = cms.string("PFAlgoDummy"),
    pfAlgoParameters = cms.PSet(
        nCalo = cms.uint32(18), 
        nMu = cms.uint32(4), # unused
        debug = cms.untracked.bool(False)
    ),
    puAlgo = cms.string("LinearizedPuppi"),
    puAlgoParameters = cms.PSet(
        nTrack = cms.uint32(0), # unused
        nIn = cms.uint32(18), 
        nOut = cms.uint32(18), 
        nVtx = cms.uint32(1),
        nFinalSort = cms.uint32(10), # to be tuned
        finalSortAlgo = cms.string("Insertion"),
        dZ     = cms.double(1.33),
        dr     = cms.double(0.3),
        drMin  = cms.double(0.1),
        ptMax  = cms.double(100.),
        absEtaCuts         = cms.vdouble(   ), # just one bin
        ptCut             = cms.vdouble( 10.0  ),
        ptSlopes           = cms.vdouble(  0.25 ),
        ptSlopesPhoton    = cms.vdouble(  0.25 ),
        ptZeros            = cms.vdouble( 14.0  ),
        ptZerosPhoton     = cms.vdouble( 14.0  ),
        alphaSlopes        = cms.vdouble(  0.6  ),
        alphaZeros         = cms.vdouble(  9.0  ),
        alphaCrop         = cms.vdouble(   4   ),
        priors             = cms.vdouble(  6.0  ),
        priorsPhoton      = cms.vdouble(  6.0  ),
        debug = cms.untracked.bool(False),
        useAssociationNetwork = cms.bool(True), #Enable Association Network
        associationThreshold = cms.double(0.1), #Association Network threshold for PV tracks
        associationGraph = cms.FileInPath("L1Trigger/L1TTrackMatch/data/NNVtx_AssociationModelGraph.pb"),
        associationNetworkZ0binning = l1tVertexProducer.VertexReconstruction.FH_HistogramParameters, #Z0 binning used for setting the input feature digitisation
        associationNetworkEtaBounds = cms.vdouble(0.0, 0.01984126984126984, 0.03968253968253968, 0.05952380952380952, 0.07936507936507936, 0.0992063492063492, 0.11904761904761904, 0.1388888888888889, 0.15873015873015872, 0.17857142857142855, 0.1984126984126984, 0.21825396825396826, 0.23809523809523808, 0.2579365079365079, 0.2777777777777778, 0.2976190476190476, 0.31746031746031744, 0.33730158730158727, 0.3571428571428571, 0.376984126984127, 0.3968253968253968, 0.41666666666666663, 0.4365079365079365, 0.45634920634920634, 0.47619047619047616, 0.496031746031746, 0.5158730158730158, 0.5357142857142857, 0.5555555555555556, 0.5753968253968254, 0.5952380952380952, 0.615079365079365, 0.6349206349206349, 0.6547619047619048, 0.6746031746031745, 0.6944444444444444, 0.7142857142857142, 0.7341269841269841, 0.753968253968254, 0.7738095238095237, 0.7936507936507936, 0.8134920634920635, 0.8333333333333333, 0.8531746031746031, 0.873015873015873, 0.8928571428571428, 0.9126984126984127, 0.9325396825396824, 0.9523809523809523, 0.9722222222222222, 0.992063492063492, 1.0119047619047619, 1.0317460317460316, 1.0515873015873016, 1.0714285714285714, 1.0912698412698412, 1.1111111111111112, 1.130952380952381, 1.1507936507936507, 1.1706349206349205, 1.1904761904761905, 1.2103174603174602, 1.23015873015873, 1.25, 1.2698412698412698, 1.2896825396825395, 1.3095238095238095, 1.3293650793650793, 1.349206349206349, 1.369047619047619, 1.3888888888888888, 1.4087301587301586, 1.4285714285714284, 1.4484126984126984, 1.4682539682539681, 1.488095238095238, 1.507936507936508, 1.5277777777777777, 1.5476190476190474, 1.5674603174603174, 1.5873015873015872, 1.607142857142857, 1.626984126984127, 1.6468253968253967, 1.6666666666666665, 1.6865079365079365, 1.7063492063492063, 1.726190476190476, 1.746031746031746, 1.7658730158730158, 1.7857142857142856, 1.8055555555555554, 1.8253968253968254, 1.8452380952380951, 1.865079365079365, 1.8849206349206349, 1.9047619047619047, 1.9246031746031744, 1.9444444444444444, 1.9642857142857142, 1.984126984126984, 2.003968253968254, 2.0238095238095237, 2.0436507936507935, 2.0634920634920633, 2.083333333333333, 2.1031746031746033, 2.123015873015873, 2.142857142857143, 2.1626984126984126, 2.1825396825396823, 2.202380952380952, 2.2222222222222223, 2.242063492063492, 2.261904761904762, 2.2817460317460316, 2.3015873015873014, 2.321428571428571, 2.341269841269841, 2.361111111111111, 2.380952380952381, 2.4007936507936507, 2.4206349206349205, 2.4404761904761902, 2.46031746031746, 2.4801587301587302, 2.5), #Eta bounds used to set z0 resolution input feature
        associationNetworkZ0ResBins = cms.vdouble(127.0, 126.0, 126.0, 126.0, 125.0, 124.0, 123.0, 122.0, 120.0, 119.0, 117.0, 115.0, 114.0, 112.0, 110.0, 107.0, 105.0, 103.0, 101.0, 98.0, 96.0, 94.0, 91.0, 89.0, 87.0, 85.0, 82.0, 80.0, 78.0, 76.0, 74.0, 72.0, 70.0, 68.0, 66.0, 64.0, 62.0, 61.0, 59.0, 57.0, 56.0, 54.0, 53.0, 51.0, 50.0, 48.0, 47.0, 46.0, 45.0, 43.0, 42.0, 41.0, 40.0, 39.0, 38.0, 37.0, 36.0, 35.0, 34.0, 33.0, 33.0, 32.0, 31.0, 30.0, 30.0, 29.0, 28.0, 28.0, 27.0, 26.0, 26.0, 25.0, 24.0, 24.0, 23.0, 23.0, 22.0, 22.0, 21.0, 21.0, 21.0, 20.0, 20.0, 19.0, 19.0, 18.0, 18.0, 18.0, 17.0, 17.0, 17.0, 16.0, 16.0, 16.0, 15.0, 15.0, 15.0, 15.0, 14.0, 14.0, 14.0, 14.0, 13.0, 13.0, 13.0, 13.0, 12.0, 12.0, 12.0, 12.0, 12.0, 11.0, 11.0, 11.0, 11.0, 11.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 9.0, 9.0, 9.0, 9.0, 9.0, 0.0), #z0 resolution input feature bins
    ),
    tkEgAlgoParameters = tkEgAlgoParameters.clone(
        nTRACK = 5,           # to be defined
        nTRACK_EGIN = 5,          # to be defined
        nEMCALO_EGIN = 5,  # to be defined
        nEM_EGOUT = 5,        # to be defined
        doBremRecovery = True,
        writeEGSta = True),
    tkEgSorterAlgo = cms.string("Endcap"),
    tkEgSorterParameters = tkEgSorterParameters.clone(),
    caloSectors = cms.VPSet(
        cms.PSet( 
            etaBoundaries = cms.vdouble(-5.5, -3.0),
            phiSlices     = cms.uint32(9),
        ),
        cms.PSet( 
            etaBoundaries = cms.vdouble(+3.0, +5.5),
            phiSlices     = cms.uint32(9),
        )
    ),
    regions = cms.VPSet(
        cms.PSet( 
            etaBoundaries = cms.vdouble(-5.5, -3.0),
            phiSlices     = cms.uint32(9),
        ),
        cms.PSet( 
            etaBoundaries = cms.vdouble(+3.0, +5.5),
            phiSlices     = cms.uint32(9),
        )
    ),
    boards = cms.VPSet(),
)


l1tLayer1 = cms.EDProducer("L1TPFCandMultiMerger",
    pfProducers = cms.VInputTag(
        cms.InputTag("l1tLayer1Barrel"),
        cms.InputTag("l1tLayer1HGCal"),
        cms.InputTag("l1tLayer1HGCalNoTK"),
        cms.InputTag("l1tLayer1HF")
    ),
)


l1tLayer1Extended = l1tLayer1.clone(
    pfProducers = [ ("l1tLayer1BarrelExtended"), ("l1tLayer1HGCalExtended"), 
        ("l1tLayer1HGCalNoTK"),("l1tLayer1HF")]
)

l1tLayer1EG = cms.EDProducer(
    "L1TEGMultiMerger",
    tkElectrons = cms.VPSet(
        cms.PSet(
            instance = cms.string("L1TkEleEE"),
            pfProducers = cms.VInputTag(
                cms.InputTag("l1tLayer1HGCal", 'L1TkEle')
            )
        ),
        cms.PSet(
            instance = cms.string("L1TkEleEB"),
            pfProducers = cms.VInputTag(
                cms.InputTag("l1tLayer1Barrel", 'L1TkEle')
            )
        )
    ),
    tkEms = cms.VPSet(
        cms.PSet(
            instance = cms.string("L1TkEmEE"),
            pfProducers = cms.VInputTag(
                cms.InputTag("l1tLayer1HGCal", 'L1TkEm'),
                cms.InputTag("l1tLayer1HGCalNoTK", 'L1TkEm')
            )
        ),
        cms.PSet(
            instance = cms.string("L1TkEmEB"),
            pfProducers = cms.VInputTag(
                cms.InputTag("l1tLayer1Barrel", 'L1TkEm')
            )
        )
    ),
    tkEgs = cms.VPSet(
        cms.PSet(
            instance = cms.string("L1EgEE"),
            pfProducers = cms.VInputTag(
                cms.InputTag("l1tLayer1HGCal", 'L1Eg'),
                cms.InputTag("l1tLayer1HGCalNoTK", 'L1Eg')
            )
        )    
    )
)

l1tLayer1EGElliptic = cms.EDProducer(
    "L1TEGMultiMerger",
    tkElectrons = cms.VPSet(
        cms.PSet(
            instance = cms.string("L1TkEleEE"),
            pfProducers = cms.VInputTag(
                cms.InputTag("l1tLayer1HGCalElliptic", 'L1TkEle')
            )
        ),
        cms.PSet(
            instance = cms.string("L1TkEleEB"),
            pfProducers = cms.VInputTag(
                cms.InputTag("l1tLayer1Barrel", 'L1TkEle')
            )
        )
    ),
    tkEms = cms.VPSet(
        cms.PSet(
            instance = cms.string("L1TkEmEE"),
            pfProducers = cms.VInputTag(
                cms.InputTag("l1tLayer1HGCalElliptic", 'L1TkEm'),
                cms.InputTag("l1tLayer1HGCalNoTK", 'L1TkEm')
            )
        ),
        cms.PSet(
            instance = cms.string("L1TkEmEB"),
            pfProducers = cms.VInputTag(
                cms.InputTag("l1tLayer1Barrel", 'L1TkEm')
            )
        )
    ),
    tkEgs = cms.VPSet(
        cms.PSet(
            instance = cms.string("L1EgEE"),
            pfProducers = cms.VInputTag(
                cms.InputTag("l1tLayer1HGCalElliptic", 'L1Eg'),
                cms.InputTag("l1tLayer1HGCalNoTK", 'L1Eg')
            )
        )    
    )
)



L1TLayer1TaskInputsTask = cms.Task(
    l1tPFClustersFromL1EGClusters,
    l1tPFClustersFromCombinedCaloHCal,
    l1tPFClustersFromCombinedCaloHF,
    l1tPFClustersFromHGC3DClusters,
    l1tPFTracksFromL1Tracks,
    l1tPFTracksFromL1TracksExtended
)

L1TLayer1Task = cms.Task(
     l1tLayer1Barrel,
     l1tLayer1BarrelExtended,
     l1tLayer1HGCal,
     l1tLayer1HGCalExtended,
     l1tLayer1HGCalNoTK,
     l1tLayer1HF,
     l1tLayer1,
     l1tLayer1Extended,
     l1tLayer1HGCalElliptic,
     l1tLayer1EG,
     l1tLayer1EGElliptic
)
