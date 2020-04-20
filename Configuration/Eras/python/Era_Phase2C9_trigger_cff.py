import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
from Configuration.Eras.Modifier_phase2_trigger_cff import phase2_trigger

Phase2C9_trigger = cms.ModifierChain(Phase2C9, phase2_trigger)

