import FWCore.ParameterSet.Config as cms
from math import pi

layer2ClusteringFw_Params = cms.PSet(
    # Latencies
    # histogramOffset=cms.uint32(216),
    clusterizerOffset=cms.uint32(276),
    stepLatencies=cms.vuint32(  0,      # Input 
                                1,      # Dist0 
                                1,      # Dist1 
                                6,      # Dist2 
                                0,      # Dist3 
                                1,      # Dist4 
                                7,      # Dist5 
                                1,      # TcToHc
                                216 + 1,# Hist, histogramOffset + 1
                                5,      # Smearing1D
                                3,      # NormArea  
                                6,      # Smearing2D
                                0,      # Maxima1D  
                                2,      # Maxima2D  
                                6,      # CalcAverage
                                0,      # Clusterizer
                                18,     # TriggerCellToCluster
                                0 ),    # ClusterSum
    cClocks=cms.uint32(250),
    cInputs=cms.uint32(72),
    cInputs2=cms.uint32(75),
    cInt=cms.uint32(90),

    minClusterPtOut=cms.double(0.5),

    # Digitzation parameters
    digiParams = cms.PSet(
        rOverZRange = cms.double(0.7),
        rOverZNValues = cms.double(4096), # 12 bits
        phiRange = cms.double(pi),
        phiNValues = cms.double(1944), # 12 bits, but range over 2pi is 3888
        ptDigiFactor = cms.double(10000),
    ),

    inputLinkParams = cms.PSet(
        maxClustersPerLink = cms.uint32(210),
        nInputLinks = cms.uint32(72)
    ),

    # Trigger cell distribution parameters
    tcDistParams = cms.PSet(
        n60Sectors = cms.uint32(3),
        nCoarsePhiRegionsDist1 = cms.uint32(6),
        nDistServers1 = cms.uint32(15),
        distServer1_nIn = cms.uint32(5),
        distServer1_nOut = cms.uint32(6),
        distServer1_nInterleave = cms.uint32(4),
        nCoarsePhiRegionsDist2 = cms.uint32(6),
        nDistServers2 = cms.uint32(18),
        distServer2_nIn = cms.uint32(5),
        distServer2_nOut = cms.uint32(6),
        distServer2_nInterleave = cms.uint32(4),
    ),

    # Histogram parameters
    cColumns=cms.uint32(108),
    cRows=cms.uint32(44),
    rOverZHistOffset = cms.uint32(440), # offset of first r/z bin in number of LSB
    rOverZBinSize = cms.uint32(64), # in number of LSB

    # Smearing parameters
    smearingParams = cms.PSet(
        maxBinsSmearing1D=cms.uint32(6),
        nBitsAreaNormLUT=cms.uint32(17)
    ),

    # Threshold maxima parameters
    # Threshold for given histogram row (r/z bin) paramaterized as a+b*bin+c*bin^2
    thresholdMaximaParams = cms.PSet(
        a=cms.uint32(18000),
        b=cms.uint32(800),
        c=cms.int32(-20)
    ),

    # Clusterizer parameters
    clusterizerParams = cms.PSet(
        nBinsCosLUT = cms.uint32(77), # Number of entries in 1-cos(delta phi) LUT when clustering to seeds
        nBitsCosLUT = cms.uint32(18), # Number of bits used in 1-cos(delta phi) LUT
        nFifos=cms.uint32(18),
        nColumnsPerFifo=cms.uint32(4),
        firstSeedBin=cms.uint32(18), # First bin considered as a seed for a cluster
        nColumnsFifoVeto=cms.uint32(6), # Minimum distance between seeds in neighbouring fifos)
        deltaR2Cut=cms.uint32(20000), # Delta R^2 threshold, in digitised units
        nColumnsForClustering=cms.uint32(3), # Number of columns either side of column containing seed to stream for clustering
        nRowsForClustering=cms.uint32(2), # Number of rows either side of row containing seed to consider during clustering
        clusterizerMagicTime=cms.uint32(434),
    ),
    depths = cms.vuint32(0 , # No zero layer
          0 , 30 , 59 , 89 , 118 , 148 , 178 , 208 , 237 , 267 , 297 , 327 , 356 , 386 , # CE-E
          415 , 445 , 475 , 505 , 534 , 564 , 594 , 624 , 653 , 683 , 712 , 742 , 772 , 802 , # CE-E
          911 , 1020 , 1129 , 1238 , 1347 , 1456 , 1565 , 1674 , # CE-FH
          1783 , 1892 , 2001 , 2110 , 2281 , 2452 , 2623 , 2794 , 2965 , 3136 , 3307 , 3478 , 3649 , 3820 , # CE-BH
          0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ),
    layerWeights_E = cms.vuint32(0 , # No zero layer
                 1 , 1 , 1 , # CE-E (early)
                 1 , 1 , 1 , 1 , 1 , # CE-E (core)
                 1 , 1 , 1 , 1 , 1 , 1 , # CE-E (back)
                 1 , 1 , 1 , 1 , # CE-H (early)
                 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1), # CE-H (back)
    layerWeights_E_EM = cms.vuint32(0 , # No zero layer
                    252969 , 254280 , 255590 , # CE-E (early)
                    256901 , 258212 , 259523 , 260833 , 262144 , # CE-E (core)
                    263455 , 264765 , 266076 , 267387 , 268698 , 270008 , # CE-E (back)
                    0 , 0 , 0 , 0 , # CE-H (early)
                    0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ), # CE-H (back)
    layerWeights_E_EM_core = cms.vuint32(0 , # No zero layer
                         0 , 0 , 0 , # CE-E (early)
                         256901 , 258212 , 259523 , 260833 , 262144 , # CE-E (core)
                         0 , 0 , 0 , 0 , 0 , 0 , # CE-E (back)
                         0 , 0 , 0 , 0 , # CE-H (early)
                         0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ), # CE-H (back)
    layerWeights_E_H_early = cms.vuint32( 0 , # No zero layer
                         0 , 0 , 0 , # CE-E (early)
                         0 , 0 , 0 , 0 , 0 , # CE-E (core)
                         0 , 0 , 0 , 0 , 0 , 0 , # CE-E (back)
                         1 , 1 , 1 , 1 , # CE-H (early)
                         0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0), # CE-H 
    correction=cms.uint32(131071), # 0b011111111111111111
    saturation=cms.uint32((2 ** 19) - 1) # 524287
)