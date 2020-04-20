import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C10_cff import Phase2C10
from Configuration.Eras.Modifier_phase2_trigger_cff import phase2_trigger

Phase2C10_trigger = cms.ModifierChain(Phase2C10, phase2_trigger)

