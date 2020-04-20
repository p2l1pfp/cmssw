import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('REPR',eras.Phase2C4_trigger)
#process = cms.Process('REPR',eras.Phase2C4_timing_layer_bar)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D35Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D35_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#'/store/mc/PhaseIIMTDTDRAutumn18DR/BsToPhiPhi_SoftQCDnonD_TuneCP5_14TeV_Pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v2/40000/9F76D756-EE13-5A48-B421-3E63846CB28F.root',
#'/store/mc/PhaseIIMTDTDRAutumn18DR/BsToPhiPhi_SoftQCDnonD_TuneCP5_14TeV_Pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v2/40000/A8265D4F-1B07-4A49-9558-2776F06DF573.root',
#'/store/mc/PhaseIIMTDTDRAutumn18DR/BsToPhiPhi_SoftQCDnonD_TuneCP5_14TeV_Pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v2/40000/B76B7085-A5A1-8B4B-946C-CCBD50C8A190.root',
#'/store/mc/PhaseIIMTDTDRAutumn18DR/BsToPhiPhi_SoftQCDnonD_TuneCP5_14TeV_Pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v2/40000/6CCEC865-6698-6B48-A20B-0C73896D680E.root',
#'/store/mc/PhaseIIMTDTDRAutumn18DR/BsToPhiPhi_SoftQCDnonD_TuneCP5_14TeV_Pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v2/40000/752602DE-CE34-5445-8445-65C1673D84DD.root',
#'/store/mc/PhaseIIMTDTDRAutumn18DR/BsToPhiPhi_SoftQCDnonD_TuneCP5_14TeV_Pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v2/40000/7540BEC9-B7C7-044D-92DF-94669EC61E68.root'
#
#'/store/mc/PhaseIIMTDTDRAutumn18DR/GluGluHToTauTau_M125_14TeV_powheg_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v1/80000/0A733B44-72DF-3A4A-A6DA-BB16EC3FDDDB.root'

'/store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v2/90000/AA56B328-0E46-E64D-93E8-4ADDF46ACEFB.root',
#'/store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v2/90000/AFD1170A-6DAA-F646-B953-9F3013BCF086.root',
#'/store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v2/90001/6F79D86E-AB09-0744-8E8E-75DE6E1022B5.root',
#'/store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v2/90001/CEF8659E-6112-5B40-9440-71F08BD2A983.root',
#'/store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v2/90001/8DFE5607-CFD8-8649-8647-F1FC67134256.root',
#'/store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v2/90001/BE6D6DDE-0274-FE4B-A6FC-E6491315DE01.root',
#'/store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v2/90001/FD13FA0E-6FB3-594D-8280-C135461881AA.root',
#'/store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v2/90001/5B441467-3183-4A42-AE41-26514EFFECAB.root',
#'/store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v2/90003/CAC9569F-AE1B-6A49-AC5F-805CDC57D529.root',
#'/store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v2/90003/409D7816-E20E-9142-9955-CCC52BAF2508.root',
#'/store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v2/90003/1468130E-7991-0E46-91A4-0B085996407F.root',
#'/store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v2/90003/8D29E3C3-69E9-F04C-A73E-5392C8522BE6.root',
#'/store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v2/90000/8DBF6AC6-0419-7448-BB32-7EB874F5A9CD.root',
#'/store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v2/90002/6D114E1A-5A35-584B-A372-DBA0F2F79DC4.root',
#'/store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v2/90003/1FE0664B-4508-1047-90D1-4F7D6819B4CC.root',
#'/store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v2/90003/1C13ED82-111C-014C-96DF-3BF6A55CBBD8.root'
),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('repr nevts:2'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
#    dataset = cms.untracked.PSet(
#        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW'),
#        filterName = cms.untracked.string('')
#    ),
    fileName = cms.untracked.string('file:step2_2ev_reprocess_slim.root'),
#    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
 outputCommands = cms.untracked.vstring('drop *', 'keep *_*_*_REPR'),
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition


# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '103X_upgrade2023_realistic_v2', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.GlobalTag = GlobalTag(process.GlobalTag, '103X_upgrade2023_realistic_v2', '') 

process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.load('CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi')

# Path and EndPath definitions
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

process.load("L1Trigger.L1TNtuples.l1PhaseIITreeProducer_cfi")


process.TFileService = cms.Service("TFileService",
    fileName = cms.string('L1NtuplePhaseII_MTD.root')
)


# Schedule definition
process.schedule = cms.Schedule(process.L1simulation_step,#process.extraCollectionsMenuTree,
            process.runmenutree,process.endjob_step)#,process.FEVTDEBUGHLToutput_step)

# Schedule definition
#process.schedule = cms.Schedule(process.L1simulation_step,process.endjob_step,process.FEVTDEBUGHLToutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)


# Customisation from command line

from L1Trigger.Configuration.customiseUtils import L1TrackTriggerTracklet
process = L1TrackTriggerTracklet(process)

process.MessageLogger.cout = cms.untracked.PSet(
    threshold = cms.untracked.string('ERROR'),
    default = cms.untracked.PSet(
        limit = cms.untracked.int32(10000)
    )
)

from L1Trigger.L1TMuonEndCap.customise_Phase2 import customise as customise_Phase2
process = customise_Phase2(process)


# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
