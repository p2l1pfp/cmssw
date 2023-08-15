import FWCore.ParameterSet.Config as cms 

from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
process = cms.Process('DIGI',Phase2C9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedHLLHC14TeV_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

# Input source
process.source = cms.Source("PoolSource",
       fileNames = cms.untracked.vstring(
        # '/store/mc/Phase2HLTTDRWinter20DIGI/SingleElectron_PT2to200/GEN-SIM-DIGI-RAW/PU200_110X_mcRun4_realistic_v3_ext2-v2/40000/00582F93-5A2A-5847-8162-D81EE503500F.root'
        'file:/hdfs/user/ec6821/HGC/LocalInputs/TTbar_200PU_V11_HLTTDR.root'
        # 'root://cms-xrd-global.cern.ch//store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/DoubleElectron_FlatPt-1To100/GEN-SIM-DIGI-RAW-MINIAOD/PU200_111X_mcRun4_realistic_T15_v1-v2/280000/003B8BCB-93B0-4040-854A-04C77E4BD066.root'
        ),
       inputCommands=cms.untracked.vstring(
           'keep *',
           'drop l1tEMTFHit2016Extras_simEmtfDigis_CSC_HLT',
           'drop l1tEMTFHit2016Extras_simEmtfDigis_RPC_HLT',
           'drop l1tEMTFHit2016s_simEmtfDigis__HLT',
           'drop l1tEMTFTrack2016Extras_simEmtfDigis__HLT',
           'drop l1tEMTFTrack2016s_simEmtfDigis__HLT',
           'drop FTLClusteredmNewDetSetVector_mtdClusters_FTLBarrel_RECO',
           'drop FTLClusteredmNewDetSetVector_mtdClusters_FTLEndcap_RECO',
           'drop MTDTrackingRecHitedmNewDetSetVector_mtdTrackingRecHits__RECO',
           'drop BTLDetIdBTLSampleFTLDataFrameTsSorted_mix_FTLBarrel_HLT',
           'drop ETLDetIdETLSampleFTLDataFrameTsSorted_mix_FTLEndcap_HLT',
           ),
        # skipEvents=cms.untracked.uint32(37)
       )

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('SingleElectronPt10_cfi nevts:10'),
    name = cms.untracked.string('Applications')
)

# Output definition
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("ntuple.root")
    )

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T15', '')

# load HGCAL TPG simulation
process.load('L1Trigger.L1THGCal.hgcalTriggerPrimitives_cff')

# Use new processors and standalone algorithms
from L1Trigger.L1THGCal.customNewProcessors import custom_clustering_standalone, custom_tower_standalone
process = custom_clustering_standalone(process)
process = custom_tower_standalone(process)

process.hgcl1tpg_step = cms.Path(process.L1THGCalTriggerPrimitives)

# Change to custom geometry
from L1Trigger.L1THGCal.customTriggerGeometry import custom_geometry_V11_Imp3
process = custom_geometry_V11_Imp3(process)

# load ntuplizer
process.load('L1Trigger.L1THGCalUtilities.hgcalTriggerNtuples_cff')
from L1Trigger.L1THGCalUtilities.customNtuples import custom_ntuples_standalone_clustering, custom_ntuples_standalone_tower
process = custom_ntuples_standalone_clustering(process)
process = custom_ntuples_standalone_tower(process)
process.ntuple_step = cms.Path(process.L1THGCalTriggerNtuples)

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    fileName = cms.untracked.string('output_withEmuHGCalClusters.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW')
    ),
)
process.FEVTDEBUGoutput.outputCommands.append( 'drop *')
process.FEVTDEBUGoutput.outputCommands.append( 'keep *_hgcalBackEndLayer2Producer_*_*')

process.FEVTDEBUGoutput.outputCommands.append( 'drop TrackingParticles_mix_MergedTrackTruth_*')
process.FEVTDEBUGoutput.outputCommands.append( 'drop PixelDigiSimLinkedmDetSetVector_simSiPixelDigis_Tracker_*')
process.FEVTDEBUGoutput.outputCommands.append( 'drop PixelDigiSimLinkedmDetSetVector_simSiPixelDigis_Pixel_*')
process.FEVTDEBUGoutput.outputCommands.append( 'drop SimClusters_mix_MergedCaloTruth_*')
process.FEVTDEBUGoutput.outputCommands.append( 'drop TrackingVertexs_mix_MergedTrackTruth_*')
process.FEVTDEBUGoutput.outputCommands.append( 'drop PixelDigiedmDetSetVector_simSiPixelDigis_Pixel_*')
process.FEVTDEBUGoutput.outputCommands.append( 'drop PCaloHits_g4SimHits_EcalHitsEB_*')
process.endjob_step = cms.EndPath(process.FEVTDEBUGoutput)

# Schedule definition
process.schedule = cms.Schedule(process.hgcl1tpg_step, process.ntuple_step, process.endjob_step )

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

