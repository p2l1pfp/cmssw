import FWCore.ParameterSet.Config as cms


tkEgSorterParameters = cms.PSet(
    splitEGObjectsToBoards = cms.untracked.bool(False),
    nObjToSort = cms.uint32(6),
    nObjSorted = cms.uint32(54),
    nBoards = cms.uint32(5)
)
