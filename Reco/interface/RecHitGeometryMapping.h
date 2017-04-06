#ifndef HGCAL_RECO_RECHITGEOMETRYMAPPING_H
#define HGCAL_RECO_RECHITGEOMETRYMAPPING_H
// user include files
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"

using namespace std;
//
// class declaration
//



class RecHitGeometryMapping
{

public:
         RecHitGeometryMapping(const edm::ParameterSet& cfg);
	~RecHitGeometryMapping();
   

	int getTBLayer(int layer);
	int getCMSSWLayer(int TBlayer);

	unsigned int getTBID_fromCMSSWID(unsigned int rawiddi);
	unsigned int getCMSSWID_fromTBID(unsigned int rawiddiTB);


private:
	int _layers_config;
	std::string _treeName;
	std::string _treeNameTB;

	std::vector<int> CMSSW_cellLayer; 

	std::map<unsigned int, unsigned int> CMSSWIDk_TBIDv;
};

#endif
