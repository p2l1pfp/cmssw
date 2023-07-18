import argparse
import sys

# example: cmsRun L1Trigger/Phase2L1ParticleFlow/test/make_l1ct_patternFiles_cfg.py -- --patternFilesOFF
# example: cmsRun L1Trigger/Phase2L1ParticleFlow/test/make_l1ct_patternFiles_cfg.py -- --dumpFilesOFF --serenity

parser = argparse.ArgumentParser(prog=sys.argv[0], description='Optional parameters')

parser.add_argument("--dumpFilesOFF", help="switch on dump file production", action="store_true", default=False)
parser.add_argument("--patternFilesOFF", help="switch on Layer-1 pattern file production", action="store_true", default=False)
parser.add_argument("--serenity", help="use Serenity settigns as default everwhere, i.e. also for barrel", action="store_true", default=False)

argv = sys.argv[:]
if '--' in argv:
    argv.remove("--")
args, unknown = parser.parse_known_args(argv)

if args.dumpFilesOFF:
    print(f'Switching off dump file creation: dumpFilesOFF is {args.dumpFilesOFF}')
if args.patternFilesOFF:
    print(f'Switching off pattern file creation: patternFilesOFF is {args.patternFilesOFF}')


import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process("RESP", eras.Phase2C17I13M9)

process.load('Configuration.StandardSequences.Services_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True), allowUnscheduled = cms.untracked.bool(False) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(480))
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:inputs110X.root'),
    inputCommands = cms.untracked.vstring("keep *", 
            "drop l1tPFClusters_*_*_*",
            "drop l1tPFTracks_*_*_*",
            "drop l1tPFCandidates_*_*_*",
            "drop l1tTkPrimaryVertexs_*_*_*")
)

process.load('Configuration.Geometry.GeometryExtended2026D88Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D88_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff') # needed to read HCal TPs
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '125X_mcRun4_realistic_v2', '')

process.load('L1Trigger.Phase2L1ParticleFlow.l1ctLayer1_cff')
process.load('L1Trigger.Phase2L1ParticleFlow.l1ctLayer2EG_cff')
process.load('L1Trigger.L1TTrackMatch.l1tGTTInputProducer_cfi')
process.load('L1Trigger.VertexFinder.l1tVertexProducer_cfi')
process.l1tVertexFinderEmulator = process.l1tVertexProducer.clone()
process.l1tVertexFinderEmulator.VertexReconstruction.Algorithm = "fastHistoEmulation"
process.l1tVertexFinderEmulator.l1TracksInputTag = cms.InputTag("l1tGTTInputProducer", "Level1TTTracksConverted")
from L1Trigger.Phase2L1GMT.gmt_cfi import l1tStandaloneMuons
process.l1tSAMuonsGmt = l1tStandaloneMuons.clone()

from L1Trigger.Phase2L1ParticleFlow.l1tSeedConePFJetEmulatorProducer_cfi import l1tSeedConePFJetEmulatorProducer
from L1Trigger.Phase2L1ParticleFlow.l1tDeregionizerProducer_cfi import l1tDeregionizerProducer
from L1Trigger.Phase2L1ParticleFlow.l1tJetFileWriter_cfi import l1tSeededConeJetFileWriter
process.l1tLayer2Deregionizer = l1tDeregionizerProducer.clone()
process.l1tLayer2SeedConeJetsCorrected = l1tSeedConePFJetEmulatorProducer.clone(L1PFObjects = cms.InputTag('l1tLayer2Deregionizer', 'Puppi'),
                                                                                doCorrections = cms.bool(True),
                                                                                correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"),
                                                                                correctorDir = cms.string('L1PuppiSC4EmuJets'))

process.runPF = cms.Path( 
        process.l1tSAMuonsGmt +
        process.l1tGTTInputProducer +
        process.l1tVertexFinderEmulator +
        process.l1tLayer1Barrel +
        process.l1tLayer1HGCal +
        process.l1tLayer1HGCalElliptic +
        process.l1tLayer1HGCalNoTK +
        process.l1tLayer1HF +
        process.l1tLayer1 +
        process.l1tLayer2Deregionizer +
        process.l1tLayer2SeedConeJetsCorrected +
        # process.l1tLayer2SeedConeJetWriter +
        process.l1tLayer2EG
    )
process.runPF.associate(process.L1TLayer1TaskInputsTask)

#####################################################################################################################
##  HGCal TPs
process.load('L1Trigger.L1THGCalUtilities.stage2FileWriter_cff')
process.Stage2FileWriter.fileExtension = "txt.gz"
from FWCore.Modules.preScaler_cfi import preScaler
for tmSlice, psOffset in (0,1), (6,2), (12,0):
    setattr(process, f"preTM{tmSlice}", preScaler.clone(prescaleFactor = 3, prescaleOffset = psOffset))
    setattr(process, f"writerTM{tmSlice}", process.Stage2FileWriter.clone(tmIndex = tmSlice))
    setattr(process, f"Write_TM{tmSlice}", cms.Path(getattr(process, f"preTM{tmSlice}")+getattr(process, f"writerTM{tmSlice}")))

# make also TM18 version.
process.Stage2FileWriter.tmIndex = 18 # a bit of a hack
process.Write_TM18 = cms.Path(process.Stage2FileWriter)

#####################################################################################################################
## Layer 1 
from L1Trigger.Phase2L1ParticleFlow.l1ctLayer1_patternWriters_cff import *
from L1Trigger.Phase2L1ParticleFlow.l1ctLayer1_patternWriters_cff import _eventsPerFile
process.l1tLayer1Barrel.patternWriters = cms.untracked.VPSet(*barrelWriterConfigs)
process.l1tLayer1HGCal.patternWriters = cms.untracked.VPSet(*hgcalWriterConfigs)
process.l1tLayer1HGCalNoTK.patternWriters = cms.untracked.VPSet(*hgcalNoTKWriterConfigs)
process.l1tLayer1HF.patternWriters = cms.untracked.VPSet(*hfWriterConfigs)

for det in "HGCal", "HGCalNoTK":
    l1pf = getattr(process, 'l1tLayer1'+det)
    l1pf.dumpFileName = cms.untracked.string("TTbar_PU200_"+det+".dump")

#####################################################################################################################
## Layer 2 e/gamma 
process.l1tLayer2EG.writeInPattern = True
process.l1tLayer2EG.writeOutPattern = True
process.l1tLayer2EG.inPatternFile.maxLinesPerFile = _eventsPerFile*54
process.l1tLayer2EG.outPatternFile.maxLinesPerFile = _eventsPerFile*54

#####################################################################################################################
## Layer 2 seeded-cone jets 
process.l1tLayer2SeedConeJetWriter = l1tSeededConeJetFileWriter.clone(jets = "l1tLayer2SeedConeJetsCorrected")
process.runPF.insert(process.runPF.index(process.l1tLayer2SeedConeJetsCorrected)+1, process.l1tLayer2SeedConeJetWriter)
process.l1tLayer2SeedConeJetWriter.maxLinesPerFile = _eventsPerFile*54

process.source.fileNames  = [ 'file:/afs/cern.ch/work/e/ejclemen/public/ForHGCalL1T/output_withEmuHGCalClusters.root' ]
process.l1tPFClustersFromHGC3DClusters.src  = cms.InputTag("l1tHGCalBackEndLayer2Producer","HGCalBackendLayer2Processor3DClusteringSA")
process.l1tPFClustersFromCombinedCaloHF.hcalCandidates = [ cms.InputTag("l1tHGCalBackEndLayer2Producer","HGCalBackendLayer2Processor3DClusteringSA")]
process.l1tPFClustersFromHGC3DClusters.corrector = ""
process.l1tPFClustersFromHGC3DClusters.emVsPUID.wp = "-99"
process.l1tPFClustersFromHGC3DClusters.useEMInterpretation = "allKeepTot"
#process.l1tLayer1HGCal.debugHGC = cms.untracked.uint32(999999)
#process.l1tLayer1HGCal.pfAlgoParameters.debug  = True
#process.l1tLayer1HGCal.puAlgoParameters.debug  = True
