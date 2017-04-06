#ifndef HGCALRECHITPRODUCER_H
#define HGCALRECHITPRODUCER_H
/** \class Reco/plugins/HGCalRecHitProducer.h HGCalRecHitProducer HGCalRecHitProducer
	\brief

	\author Shervin Nourbakhsh
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBDigiCollections.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
//#include "HGCal/DataFormats/plugins/HGCGeometryReader.cc"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// here are defined objects containing pedestals and ADCtoGeV factors
#include "HGCal/CondObjects/interface/HGCalCondObjects.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/CondObjects/interface/HGCalTBNumberingScheme.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Reco/interface/RecHitGeometryMapping.h"


#include "TTree.h"
#include "TH2F.h"
#include <iostream>

class HGCalRecHitProducer : public edm::one::EDProducer<edm::one::SharedResources>
{

public:
	HGCalRecHitProducer(const edm::ParameterSet&);
	virtual void produce(edm::Event&, const edm::EventSetup&);

	/*
	int getCellLayer(int TBlayer){
	  //std::cout << " TBlayer = " << TBlayer << "corresponds to = " << CMSSW_cellLayer.at(TBlayer-1) << std::endl;
	  return CMSSW_cellLayer.at(TBlayer-1);
	}

	int getCellID(float lX, float lY){
	  //std::cout << " bin = " << waferMap->FindBin(lX, lY) << std::endl;
	  return waferMap->GetBinContent(waferMap->FindBin(lX, lY));
	}

	unsigned int getRawID(std::pair<int, int> layerCell){
	  return (LayerCellID[layerCell]);
	}
	*/

private:
	std::string outputCollectionName;     ///<label name of collection made by this producer
	edm::EDGetTokenT<HGCalTBRecHitCollection> _inputTBCollection;
	int _layers_config;
	//	std::string _treeName;
	//      std::string _treeNameTB;

	const int MaxNlayer = 8;
	const int sensorsize = 128;
	const int CMthreshold = 2;
	/* double commonmode[MaxNlayer]; */
	/* int cm_num[MaxNlayer]; */
	double commonmode[8];
	int cm_num[8];
	double weights2MIP = 1.; 
	double weights2GeV = 1.e-03;
	double MIP2GeV_sim = 51.91e-06;

	std::string _mapFile;
	struct {
	  HGCalElectronicsMap emap_;
        } essource_;

	std::pair<double, double> CellCentreXY;
	HGCalTBCellVertices TheCell;
	//	std::vector<int> CMSSW_cellLayer;
	std::vector<double> Weights_L;
	std::vector<double> ADCtoMIP;
 
	TH1F* h_layer;
	TH1F* h_IU;
	TH1F* h_IV;
	TH1F* h_iu;
	TH1F* h_iv;
	TH1F* h_cellType;


	TTree* TBtree;
	Int_t layerTB;
	Float_t localXTB;
	Float_t localYTB;
	UInt_t detIDTB;

	RecHitGeometryMapping* TBcmsswGeometryMap;

};



#endif
