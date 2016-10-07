/* Modify RecHitPlotter_HighGain_Correlation_CM.cc to use the RecHitCommonMode class */
//
// Author:  Menglei Sun
// Created: Fri, August 26 2016


// system include files
#include <memory>
#include <iostream>
#include "TH2Poly.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TH2F.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#include "HGCal/Reco/interface/RecHitCommonMode.h"

using namespace std;
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

static const double delta = 0.00001;//Add/subtract delta = 0.00001 to x,y of a cell centre so the TH2Poly::Fill doesnt have a problem at the edges where the centre of a half-hex cell passes through the sennsor boundary line.

class RecHitPlotter_HighGain_CM_Correction : public edm::one::EDAnalyzer<edm::one::SharedResources>
{

public:
	explicit RecHitPlotter_HighGain_CM_Correction(const edm::ParameterSet&);
	~RecHitPlotter_HighGain_CM_Correction();
	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
	virtual void beginJob() override;
	void analyze(const edm::Event& , const edm::EventSetup&) override;
	virtual void endJob() override;

	// ----------member data ---------------------------
	bool doCommonMode_CM_; 
	edm::EDGetToken HGCalTBRecHitCollection_;
	HGCalTBTopology IsCellValid;
	HGCalTBCellVertices TheCell;
        std::string mapfile_ = "HGCal/CondObjects/data/map_FNAL_SB2_Layer16.txt";
        struct {
                HGCalElectronicsMap emap_;
        } essource_;
	int sensorsize = 128;// The geometry for a 256 cell sensor hasnt been implemted yet. Need a picture to do this.
	std::vector<std::pair<double, double>> CellXY;
	std::pair<double, double> CellCentreXY;

	TH2Poly *h_RecHit_layer[128];
        TH1F* Full_Cell[MAXLAYERS];
        TH1F* Half_Cell[MAXLAYERS];
        TH1F* MB_Cell[MAXLAYERS];
        TH1F* Calib_Pads[MAXLAYERS];
        TH1F* Merged_Cell[MAXLAYERS];
        TH1F  *h_digi_layer_channel[MAXSKIROCS][64][MAXLAYERS];

        TH1F* AllCells_Ped;
        TH1F* AllCells_CM;
        TH2F* Noise_2D_Profile;   
	char name[50], title[50];
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
RecHitPlotter_HighGain_CM_Correction::RecHitPlotter_HighGain_CM_Correction(const edm::ParameterSet& iConfig)
{
	//now do what ever initialization is needed
	usesResource("TFileService");
        doCommonMode_CM_ = iConfig.getParameter<bool>("doCommonMode");
	edm::Service<TFileService> fs;
	HGCalTBRecHitCollection_ = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));


//Booking 2 "hexagonal" histograms to display the sum of Rechits and the Occupancy(Hit > 5 GeV) in 1 sensor in 1 layer. To include all layers soon. Also the 1D Rechits per cell in a sensor is booked here.
        AllCells_Ped = fs->make<TH1F>("AllCells_Ped","AllCells_Ped",500,-250,250);
        AllCells_CM = fs->make<TH1F>("AllCells_CM","AllCells_CM",500,-250,250);
        sprintf(name, "Noise_2D_Profile_Layer");
        sprintf(title, "Noise 2D Profile Layer");
        Noise_2D_Profile = fs->make<TH2F>(name,title,2048,0,2048,8000,-4000,4000);
        for(int ILayer=0;ILayer<MAXLAYERS;ILayer++){ 
             sprintf(name, "Full_Cell_Layer_%i",ILayer);
             sprintf(title, "Full Cell Layer %i",ILayer);
             Full_Cell[ILayer] = fs->make<TH1F>(name, title, 1000,-500., 500.);
             sprintf(name, "Half_Cell_Layer_%i",ILayer);
             sprintf(title, "Half Cell Layer %i",ILayer);
             Half_Cell[ILayer] = fs->make<TH1F>(name, title, 1000,-500., 500.);
             sprintf(name, "MB_Cell_Layer_%i",ILayer);
             sprintf(title, "MB Cell Layer %i",ILayer);
             MB_Cell[ILayer] = fs->make<TH1F>(name, title, 1000,-500., 500.);
             sprintf(name, "Calib_Pads_Layer_%i",ILayer);
             sprintf(title, "Calib Pads Layer %i",ILayer);
             Calib_Pads[ILayer] = fs->make<TH1F>(name, title, 1000,-500., 500.);             
             sprintf(name, "Merged_Cell_Layer_%i",ILayer);
             sprintf(title, "Merged Cell Layer %i",ILayer);
             Merged_Cell[ILayer] = fs->make<TH1F>(name, title, 1000,-500., 500.);
	        for(int ISkiroc = 1;ISkiroc<=MAXSKIROCS;ISkiroc++){
        		for(int Channel=0; Channel<64;Channel++){
				sprintf(name, "Ski_%i_Channel_%i_Layer_%i",ISkiroc,Channel,ILayer);
                  		sprintf(title, "Ski %i Channel %i Layer %i",ISkiroc,Channel,ILayer);
				h_digi_layer_channel[ISkiroc-1][Channel][ILayer] = fs->make<TH1F>(name, title, 1000,-500., 500.);
/*
        			 sprintf(name, "Ski_%i_Channel_%i_CM",ISkiroc,Channel);
       			         sprintf(title, "Ski %i Channel %i CM",ISkiroc,Channel);
                		 h_digi_layer_channel_CM[ISkiroc-1][Channel] = fs->make<TH2F>(name, title, 1000,-500., 500.,1000,-500., 500.);
*/
            			  }
          		 }
        	 }

}//contructor ends here


RecHitPlotter_HighGain_CM_Correction::~RecHitPlotter_HighGain_CM_Correction()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
RecHitPlotter_HighGain_CM_Correction::analyze(const edm::Event& event, const edm::EventSetup& setup)
{

	using namespace edm;

	edm::Handle<HGCalTBRecHitCollection> Rechits;
	event.getByToken(HGCalTBRecHitCollection_, Rechits);
        edm::Handle<HGCalTBRecHitCollection> Rechits1;
        event.getByToken(HGCalTBRecHitCollection_, Rechits1);

        RecHitCommonMode CM(Rechits);
        CM.evaluate();
        
        for(int iii = 0; iii < MAXLAYERS; iii++){
	  if(CM.getCommonModeNoise(iii, HGCalTBDetId::kCellTypeStandard) != 0) Full_Cell[iii]->Fill(CM.getCommonModeNoise(iii, HGCalTBDetId::kCellTypeStandard));
	  if(CM.getCommonModeNoise(iii, HGCalTBDetId::kCellTypeHalfCell) != 0) Half_Cell[iii]->Fill(CM.getCommonModeNoise(iii, HGCalTBDetId::kCellTypeHalfCell));
	  if(CM.getCommonModeNoise(iii, HGCalTBDetId::kCellTypeMerged, "MB") != 0) MB_Cell[iii]->Fill(CM.getCommonModeNoise(iii, HGCalTBDetId::kCellTypeMerged, "MB"));
	  if(CM.getCommonModeNoise(iii, HGCalTBDetId::kCellTypeMerged, "MergedCell") != 0) Merged_Cell[iii]->Fill(CM.getCommonModeNoise(iii, HGCalTBDetId::kCellTypeMerged, "MergedCell"));
	  if(CM.getCommonModeNoise(iii, HGCalTBDetId::kCellTypeCalibInner) != 0) Calib_Pads[iii]->Fill(CM.getCommonModeNoise(iii, HGCalTBDetId::kCellTypeCalibInner));
	  if(CM.getCommonModeNoise(iii, HGCalTBDetId::kCellTypeCalibOuter) != 0) Calib_Pads[iii]->Fill(CM.getCommonModeNoise(iii, HGCalTBDetId::kCellTypeCalibOuter));
        }        

	for(auto RecHit : *Rechits) {
          if(!IsCellValid.iu_iv_valid((RecHit.id()).layer(), (RecHit.id()).sensorIU(), (RecHit.id()).sensorIV(), (RecHit.id()).iu(), (RecHit.id()).iv(), sensorsize))  continue;
          CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots((RecHit.id()).layer(), (RecHit.id()).sensorIU(), (RecHit.id()).sensorIV(), (RecHit.id()).iu(), (RecHit.id()).iv(), sensorsize);

          uint32_t EID = essource_.emap_.detId2eid(RecHit.id());
          HGCalTBElectronicsId eid(EID);
	  if(!doCommonMode_CM_){
	    h_digi_layer_channel[eid.iskiroc()-1][eid.ichan()][(RecHit.id()).layer() -1]->Fill(RecHit.energyHigh());
	    Noise_2D_Profile->Fill((64*(eid.iskiroc()-1) + eid.ichan()),RecHit.energyHigh());  
	  }
	  else if(doCommonMode_CM_){
	    AllCells_Ped->Fill(RecHit.energyHigh());
	    h_digi_layer_channel[eid.iskiroc()-1][eid.ichan()][(RecHit.id()).layer() -1]->Fill(RecHit.energyHigh() - CM.getCommonModeNoise(RecHit.id()));
          }  
	  
          Noise_2D_Profile->Fill((64*(eid.iskiroc()-1) + eid.ichan()),RecHit.energyHigh() - CM.getCommonModeNoise((RecHit.id()).layer() -1, HGCalTBDetId::kCellTypeStandard));

	} 
                           

}//analyze method ends here


// ------------ method called once each job just before starting event loop  ------------
void
RecHitPlotter_HighGain_CM_Correction::beginJob()
{
HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(mapfile_);
   if (!io.load(fip.fullPath(), essource_.emap_)) {
     throw cms::Exception("Unable to load electronics map");
      }
}

// ------------ method called once each job just after ending the event loop  ------------
void
RecHitPlotter_HighGain_CM_Correction::endJob()
{

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RecHitPlotter_HighGain_CM_Correction::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RecHitPlotter_HighGain_CM_Correction);
