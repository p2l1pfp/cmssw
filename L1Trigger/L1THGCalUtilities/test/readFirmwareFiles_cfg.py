import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing


# PART 1 : PARSE ARGUMENTS

options = VarParsing.VarParsing ('analysis')
options.parseArguments()

inputFiles = []
for filePath in options.inputFiles:
    if filePath.endswith(".root"):
        inputFiles.append(filePath)
    else:
        inputFiles += FileUtils.loadListFromFile(filePath)


# PART 2: SETUP MAIN CMSSW PROCESS 

process = cms.Process("HGCClusterValidation")

process.load('Configuration.Geometry.GeometryExtended2026D88Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D88_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(inputFiles) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.load('L1Trigger.L1THGCalUtilities.stage2FileReader_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 3
process.Timing = cms.Service("Timing", summaryOnly = cms.untracked.bool(True))

# process.Stage2FileReader.files = cms.vstring("output_hgcAlgo.txt")

process.Stage2FileReader.files = cms.vstring(
    'output_hgcAlgo_0.txt',
    'output_hgcAlgo_1.txt',
    'output_hgcAlgo_2.txt',
    'output_hgcAlgo_3.txt',
    'output_hgcAlgo_4.txt',
    'output_hgcAlgo_5.txt'
    )


process.load('FWCore.Modules.preScaler_cfi')
process.preScaler.prescaleFactor = 3
process.preScaler.prescaleOffset = 1
process.p = cms.Path(process.preScaler*process.Stage2FileReader)
