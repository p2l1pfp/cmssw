import FWCore.ParameterSet.Config as cms
from math import pi
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing


process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

options = VarParsing.VarParsing ('analysis')
# get and parse the command line arguments

options.register('skipEvents',
                 0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Number of events to skip")
options.register('outFile',
                 'L1Ntuple.root',
                  VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file')

options.parseArguments()

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

fileList = FileUtils.loadListFromFile('snu.list')
readFiles = cms.untracked.vstring(*fileList)

process.source = process.source = cms.Source("PoolSource",
  fileNames = readFiles,
)

process.load('L1Trigger.L1CaloTrigger.Phase1L1TJets_cff')

process.load('L1Trigger.Phase2L1ParticleFlow.l1pfJetMet_cff')
process.l1PFJets.remove(process.ak4PFL1Calo)
process.l1PFJets.remove(process.ak4PFL1PF)
process.l1PFJets.remove(process.ak4PFL1CaloCorrected)
process.l1PFJets.remove(process.ak4PFL1PFCorrected)


process.load("L1Trigger.L1TNtuples.l1PhaseIPFJetTreeProducer_cfi")


process.TFileService = cms.Service("TFileService",
    fileName = cms.string('L1Ntuple.root')
)

process.p = cms.Path(process.Phase1L1TJetsSequence+process.l1PFJets)
