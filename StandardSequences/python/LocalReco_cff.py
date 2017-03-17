import FWCore.ParameterSet.Config as cms

from HGCal.Reco.hgcaltbrechitproducer_cfi import *
from HGCal.Reco.hgcalrechitproducer_cfi import *

LocalRecoSeq  = cms.Sequence(hgcaltbrechits*HGCalRecHit)


