# define basic process
import FWCore.ParameterSet.Config as cms
import os
process = cms.Process("L1Tracklet")
 

# import standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.Geometry.GeometryExtended2023D35Reco_cff') ## this needs to match the geometry used
process.load('Configuration.Geometry.GeometryExtended2023D35_cff')     ## this needs to match the geometry used
#process.load('Configuration.Geometry.GeometryExtended2023D41Reco_cff') ## this needs to match the geometry used (D41: change in ../interface/Constants.h to geomTkTDR=false)
#process.load('Configuration.Geometry.GeometryExtended2023D41_cff')     ## this needs to match the geometry used (D41: change in ../interface/Constants.h to geomTkTDR=false)

process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')


# input
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(5))
Source_Files = cms.untracked.vstring(
    "/store/mc/PhaseIIMTDTDRAutumn18DR/TTbar_14TeV_TuneCP5_Pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v1/80003/FFF35381-A845-4841-9364-923285FFCFA6.root"
#    "/store/mc/PhaseIITDRSpring19DR/TTbar_14TeV_TuneCP5_Pythia8/GEN-SIM-DIGI-RAW/PU200_106X_upgrade2023_realistic_v3_ext1-v3/60000/FFB5D0CA-208F-6040-A9BF-3F5354D0AA59.root"
    )
process.source = cms.Source("PoolSource", fileNames = Source_Files)


# L1 tracking => floating-point version
#process.load("L1Trigger.TrackFindingTracklet.L1TrackletTracks_cff")
#process.TTTracks = cms.Path(process.L1TrackletTracks)                         #run only the tracking (no MC truth associators)
#process.TTTracksWithTruth = cms.Path(process.L1TrackletTracksWithAssociators) #run the tracking AND MC truth associators)

# L1 tracking => hybrid emulation 
#process.load("L1Trigger.TrackFindingTracklet.L1TrackletEmulationTracks_cff")
#process.TTTracksEmulation = cms.Path(process.L1TrackletEmulationTracks)
#process.TTTracksEmulationWithTruth = cms.Path(process.L1TrackletEmulationTracksWithAssociators)

# L1 tracking => extended hybrid emulation (displaced)                  
#process.load("L1Trigger.TrackFindingTracklet.L1ExtendedTrackletEmulationTracks_cff")      
#process.TTTracksExtendedEmulation = cms.Path(process.L1ExtendedTrackletEmulationTracks)
#process.TTTracksExtendedEmulationWithTruth = cms.Path(process.L1ExtendedTrackletEmulationTracksWithAssociators)  

# L1 tracking => floating-point version + hybrid emulation
process.load("L1Trigger.TrackFindingTracklet.L1TrackletTracks_cff")
process.TTTracksWithTruth = cms.Path(process.L1TrackletTracksWithAssociators)
process.load("L1Trigger.TrackFindingTracklet.L1TrackletEmulationTracks_cff")
process.TTTracksEmulation = cms.Path(process.L1TrackletEmulationTracks)

# output module
process.out = cms.OutputModule( "PoolOutputModule",
                                fileName = cms.untracked.string("Tracklets.root"),
                                fastCloning = cms.untracked.bool( False ),
                                outputCommands = cms.untracked.vstring('drop *',
                                                                       'keep *_TTTrack*_Level1TTTracks_*', 
#                                                                       'keep *_TTCluster*_*_*',
#                                                                       'keep *_TTStub*_*_*'
)
)
process.FEVToutput_step = cms.EndPath(process.out)

process.schedule = cms.Schedule(process.TTTracksWithTruth,process.TTTracksEmulation,process.FEVToutput_step)

