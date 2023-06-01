import argparse
import sys

# example: cmsRun L1Trigger/Phase2L1ParticleFlow/test/make_l1ct_patternFiles_cfg.py -- --dumpFilesOFF
# example: cmsRun L1Trigger/Phase2L1ParticleFlow/test/make_l1ct_patternFiles_cfg.py -- --dumpFilesOFF

parser = argparse.ArgumentParser(prog=sys.argv[0], description='Optional parameters')

parser.add_argument("--dumpFilesOFF", help="switch on dump file production", action="store_true", default=False)
parser.add_argument("--patternFilesOFF", help="switch on Layer-1 pattern file production", action="store_true", default=False)

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

process = cms.Process("RESP", eras.Phase2C9)

process.load('Configuration.StandardSequences.Services_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True), allowUnscheduled = cms.untracked.bool(False) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(24))
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:inputs110X.root'),
    inputCommands = cms.untracked.vstring("keep *", 
            "drop l1tPFClusters_*_*_*",
            "drop l1tPFTracks_*_*_*",
            "drop l1tPFCandidates_*_*_*",
            "drop l1tTkPrimaryVertexs_*_*_*")
)

process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff') # needed to read HCal TPs
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '123X_mcRun4_realistic_v3', '')

process.load('L1Trigger.Phase2L1ParticleFlow.l1ctLayer1_cff')
process.load('L1Trigger.Phase2L1ParticleFlow.l1ctLayer2EG_cff')
process.load('L1Trigger.L1TTrackMatch.l1tGTTInputProducer_cfi')
process.load('L1Trigger.VertexFinder.l1tVertexProducer_cfi')
process.l1tVertexFinderEmulator = process.l1tVertexProducer.clone()
process.l1tVertexFinderEmulator.VertexReconstruction.Algorithm = "fastHistoEmulation"
process.l1tVertexFinderEmulator.l1TracksInputTag = cms.InputTag("l1tGTTInputProducer", "Level1TTTracksConverted")
from L1Trigger.Phase2L1GMT.gmt_cfi import l1tStandaloneMuons
process.l1tSAMuonsGmt = l1tStandaloneMuons.clone()

from L1Trigger.Phase2L1ParticleFlow.l1tSeedConePFJetProducer_cfi import l1tSeedConePFJetEmulatorProducer
from L1Trigger.Phase2L1ParticleFlow.l1tDeregionizerProducer_cfi import l1tDeregionizerProducer
from L1Trigger.Phase2L1ParticleFlow.l1tJetFileWriter_cfi import l1tSeededConeJetFileWriter
process.l1tLayer2Deregionizer = l1tDeregionizerProducer.clone()
process.l1tLayer2SeedConeJetsCorrected = l1tSeedConePFJetEmulatorProducer.clone(L1PFObject = cms.InputTag('l1tLayer2Deregionizer', 'Puppi'),
                                                                                doCorrections = cms.bool(True),
                                                                                correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"),
                                                                                correctorDir = cms.string('L1PuppiSC4EmuJets'))
process.l1tLayer2SeedConeJetWriter = l1tSeededConeJetFileWriter.clone(jets = "l1tLayer2SeedConeJetsCorrected")


from L1Trigger.Phase2L1ParticleFlow.l1ctLayer1_patternWriters_cff import hgcalWriterConfig_, hgcalNoTKWriterOutputOnlyConfig 
if not args.patternFilesOFF:
    hgcalPosVU9PWriterConfig = hgcalWriterConfig_.clone()
    for t in range(3):
        hgcalPosVU9PWriterConfig.tfTimeSlices[t].tfSectors += [ cms.PSet(tfLink = cms.int32(-1))    for i in range(9) ] # neg
        hgcalPosVU9PWriterConfig.tfTimeSlices[t].tfSectors += [ cms.PSet(tfLink = cms.int32(3*i+t+4*2)) for i in range(4) ] # pos, left quads
        hgcalPosVU9PWriterConfig.tfTimeSlices[t].tfSectors += [ cms.PSet(tfLink = cms.int32(3*i+t+4*25)) for i in range(5) ] # pos, right quads
        hgcalPosVU9PWriterConfig.hgcTimeSlices[t].hgcSectors += [ cms.PSet(hgcLinks = cms.vint32(-1,-1,-1,-1))                        for i in range(3) ] # neg
        hgcalPosVU9PWriterConfig.hgcTimeSlices[t].hgcSectors += [ cms.PSet(hgcLinks = cms.vint32(*[4*11+12*i+4*t+j for j in range(4)])) for i in range(3) ] # pos
        hgcalPosVU9PWriterConfig.gmtTimeSlices[t].gmtLink = cms.int32(4+t)
    hgcalPosVU9PWriterConfig.gttLink = 4+3
    hgcalPosVU9PWriterConfig.outputLinksPuppi = cms.vuint32(56,57,58)
    hgcalPosVU9PWriterConfig.outputLinkEgamma = cms.int32(59)
    hgcalPosVU9PWriterConfig.inputFileName = cms.string("l1HGCalPos-inputs-vu9p") 
    hgcalPosVU9PWriterConfig.outputFileName = cms.string("l1HGCalPos-outputs-vu9p")
    hgcalNegVU9PWriterConfig = hgcalWriterConfig_.clone()
    for t in range(3):
        hgcalNegVU9PWriterConfig.tfTimeSlices[t].tfSectors += [ cms.PSet(tfLink = cms.int32(3*i+t+4*2)) for i in range(4) ] # pos, left quads
        hgcalNegVU9PWriterConfig.tfTimeSlices[t].tfSectors += [ cms.PSet(tfLink = cms.int32(3*i+t+4*25)) for i in range(5) ] # pos, right quads
        hgcalNegVU9PWriterConfig.tfTimeSlices[t].tfSectors += [ cms.PSet(tfLink = cms.int32(-1))    for i in range(9) ] # neg
        hgcalNegVU9PWriterConfig.hgcTimeSlices[t].hgcSectors += [ cms.PSet(hgcLinks = cms.vint32(*[4*11+12*i+4*t+j for j in range(4)])) for i in range(3) ] # pos
        hgcalNegVU9PWriterConfig.hgcTimeSlices[t].hgcSectors += [ cms.PSet(hgcLinks = cms.vint32(-1,-1,-1,-1))                          for i in range(3) ] # neg
        hgcalNegVU9PWriterConfig.gmtTimeSlices[t].gmtLink = cms.int32(4+t)
    hgcalNegVU9PWriterConfig.gttLink = 4+3
    hgcalNegVU9PWriterConfig.outputLinksPuppi = cms.vuint32(56,57,58)
    hgcalNegVU9PWriterConfig.outputLinkEgamma = cms.int32(59)
    hgcalNegVU9PWriterConfig.inputFileName = cms.string("l1HGCalNeg-inputs-vu9p") 
    hgcalNegVU9PWriterConfig.outputFileName = cms.string("l1HGCalNeg-outputs-vu9p")

    process.l1tLayer1HGCal.patternWriters = cms.untracked.VPSet(hgcalPosVU9PWriterConfig,hgcalNegVU9PWriterConfig)
    process.l1tLayer1HGCalNoTK.patternWriters = cms.untracked.VPSet()

process.runPF = cms.Path( 
        process.l1tSAMuonsGmt +
        process.l1tGTTInputProducer +
        process.l1tVertexFinderEmulator +
        process.l1tLayer1HGCal +
        process.l1tLayer1HGCalElliptic +
        process.l1tLayer1HGCalNoTK
    )
process.runPF.associate(process.L1TLayer1TaskInputsTask)

if not args.dumpFilesOFF:
  for det in "HGCal", "HGCalElliptic", "HGCalNoTK":
        l1pf = getattr(process, 'l1tLayer1'+det)
        l1pf.dumpFileName = cms.untracked.string("TTbar_PU200_HGCTest_"+det+".dump")

### FIXES for integration tests
process.source.fileNames  = [ 'file:/afs/cern.ch/work/e/ejclemen/public/ForHGCalL1T/output_withEmuHGCalClusters.root' ]
process.l1tPFClustersFromHGC3DClusters.src  = cms.InputTag("l1tHGCalBackEndLayer2Producer","HGCalBackendLayer2Processor3DClusteringSA")
process.l1tPFClustersFromCombinedCaloHF.hcalCandidates = [ cms.InputTag("l1tHGCalBackEndLayer2Producer","HGCalBackendLayer2Processor3DClusteringSA")]
process.l1tPFClustersFromHGC3DClusters.corrector = ""
process.l1tPFClustersFromHGC3DClusters.emVsPUID.wp = "-99"
process.l1tPFClustersFromHGC3DClusters.useEMInterpretation = "allKeepTot"

if True:
    process.l1tLayer1HGCal.debugHGC = cms.untracked.uint32(2000) # dump the first 2000 clusters per job
    process.l1tLayer1HGCal.tkEgAlgoParameters.debug = cms.untracked.uint32(4)
    process.l1tLayer1HGCal.pfAlgoParameters.debug = True 

