import FWCore.ParameterSet.Config as cms

HGCalRecHit = cms.EDProducer("HGCalRecHitProducer",
                             OutputCollectionName = cms.string('HGCEERecHits'),
                             #inputTBCollection = cms.InputTag('hgcaltbrechits'),
                             inputTBCollection = cms.InputTag("hgcaltbrechits","","unpack" ),
                             layers_config = cms.int32(-1),
                             treeName = cms.string('../testHGCGeometry.root'),
                             treeNameTB = cms.string('../testTBGeometry.root'),
                             mapFile_CERN = cms.string('HGCal/CondObjects/data/map_CERN_8Layers_Sept2016.txt'),
                             mapFile_FNAL = cms.string(''), 
                             CMSSW_cellLayer_1 = cms.vint32(10, 11, 12, 13, 14, 16, 17, 19),
                             CMSSW_cellLayer_2 = cms.vint32(8, 12, 16, 19, 22, 23, 25, 28),
                             LayerWeight_8L_conf1 = cms.vdouble(33.074, 13.184, 14.17, 9.788, 9.766, 9.766, 16.339, 14.129),
                             LayerWeight_8L_conf2 = cms.vdouble(35.866, 30.864, 28.803, 23.095, 20.657, 19.804, 36.322, 27.451),
                             ADCtoMIP_CERN = cms.vdouble(17.31, 17.12, 16.37, 17.45, 17.31, 16.98, 16.45, 16.19, 17.55, 17.19, 16.99, 17.92, 15.95, 16.64, 16.79, 15.66)
                             #treeName = cms.string('HGCal/Reco/plugins/testHGCGeometry.root')
                             
                              )
