#include "HGCal/Reco/plugins/HGCalRecHitProducer.h"
//#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
//#include "DataFormats/plugins/HGCGeometryReader.cc"
#include <iostream>

#include "TTree.h"
#include "TFile.h"
#include "TH2F.h"

using namespace std;

HGCalRecHitProducer::HGCalRecHitProducer(const edm::ParameterSet& cfg): 
  outputCollectionName(cfg.getParameter<std::string>("OutputCollectionName")),
  _inputTBCollection(consumes<HGCalTBRecHitCollection>(cfg.getParameter<edm::InputTag>("inputTBCollection"))),
  _layers_config(cfg.getParameter<int>("layers_config")), 
  _treeName(cfg.getParameter<std::string>("treeName")){

  //  produces <HGCRecHitCollection>(outputCollectionName);
  //  produces<HGCeeRecHitCollection>(outputCollectionName);
  produces<edm::SortedCollection<HGCRecHit> >(outputCollectionName);
  if(_layers_config == 0){
    _mapFile = cfg.getParameter<std::string>("mapFile_FNAL");
  }
  else{
    _mapFile = cfg.getParameter<std::string>("mapFile_CERN");
    ADCtoMIP = cfg.getParameter<std::vector<double> >("ADCtoMIP_CERN");
    //fixme check order 
    if(_layers_config == 1){
      CMSSW_cellLayer = cfg.getParameter<std::vector<int> >("CMSSW_cellLayer_1");
      Weights_L = cfg.getParameter<std::vector<double> >("LayerWeight_8L_conf1");
    }
    if(_layers_config == 2) {
      CMSSW_cellLayer = cfg.getParameter<std::vector<int> >("CMSSW_cellLayer_2");
      Weights_L = cfg.getParameter<std::vector<double> >("LayerWeight_8L_conf2");
    }
  }
  //  std::cout << " HGCalRecHitProducer::HGCalRecHitProducer get File " << std::endl;

  // to load the detID 
  TFile* inF = TFile::Open(_treeName.c_str());
  TTree* tree = (TTree*)inF->Get("ana/t");
  //  std::cout << " entries = " << tree->GetEntries() << std::endl;

  Int_t   subdet, layer, wafer, cell;
  Float_t z, locX, locY;
  UInt_t detID;

  tree->SetBranchAddress("subdet", &subdet);
  tree->SetBranchAddress("layer",&layer);
  tree->SetBranchAddress("wafer",&wafer);
  tree->SetBranchAddress("cell",&cell);
  tree->SetBranchAddress("detID", &detID);
  tree->SetBranchAddress("z",&z);
  tree->SetBranchAddress("locX",&locX);
  tree->SetBranchAddress("locY",&locY);

  //  waferMap = new TH2F("waferMap", "",  200, -10, 10, 200, -10, 10);
  waferMap = new TH2F("waferMap", "", 23, -6.5, 6.5, 15, -7., 7.);
  for(int ij=0; ij<tree->GetEntries(); ++ij){
    tree->GetEntry(ij);

    //    std::cout << " layer = " << layer << " wafer = " << wafer << " z = " << z << " subdet = " << subdet << " detID = " << detID << " cellID = " << cell << " bin = " << waferMap->FindBin(locX, locY) << std::endl;
    if(wafer == 426 && z > 0 && subdet == 3 && cell < 132){    
      if(layer == 1) waferMap->SetBinContent(waferMap->FindBin(locX, locY), cell);
      std::pair<int, int> LayerCell(layer, cell);
      LayerCellID[LayerCell] = detID;
    }
  }

  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(_mapFile);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  };
  
  for(int iL=0; iL<MaxNlayer; ++iL){
      commonmode[iL] = 0.;
      cm_num[iL] = 0;
  }
}



void HGCalRecHitProducer::produce(edm::Event& event, const edm::EventSetup& iSetup)
{
  //  std::auto_ptr<HGCRecHitCollection> rechits(new HGCRecHitCollection); //auto_ptr are automatically deleted when out of scope
  //  std::auto_ptr<HGCeeRecHitCollection> rechits(new HGCeeRecHitCollection); //auto_ptr are automatically deleted when out of scope
  //std::auto_ptr<edm::SortedCollection<HGCRecHit> > rechits(new edm::SortedCollection<HGCRecHit>); //auto_ptr are automatically deleted when out of scope
  auto rechits = std::make_unique<edm::SortedCollection<HGCRecHit> >();


  edm::Handle<HGCalTBRecHitCollection> tbRecHitsHandle;
  event.getByToken(_inputTBCollection, tbRecHitsHandle);
  
  //  std::cout << " TBRH size = " << tbRecHitsHandle->size() << std::endl;

  //compute CM
  for(auto tbRH_itr = tbRecHitsHandle->begin(); tbRH_itr != tbRecHitsHandle->end(); ++tbRH_itr) {

    uint32_t EID = essource_.emap_.detId2eid(tbRH_itr->id());
    HGCalTBElectronicsId eid(EID);
    int n_skiroc = (eid.iskiroc() - 1);
    int n_layer = (tbRH_itr->id()).layer() - 1;

    int n_cell_type = (tbRH_itr->id()).cellType();
    if(n_cell_type != 0 && n_cell_type != 4) continue;

    if(tbRH_itr->energy() / ADCtoMIP[n_skiroc] <= CMthreshold) {
      commonmode[n_layer] += tbRH_itr->energy();
      cm_num[n_layer]++;
    }
  }
  for(int iL=0; iL<MaxNlayer; ++iL){
    commonmode[iL] = commonmode[iL]/cm_num[iL];
  }


  for(auto tbRH_itr = tbRecHitsHandle->begin(); tbRH_itr != tbRecHitsHandle->end(); ++tbRH_itr) {

    int n_cell_type = (tbRH_itr->id()).cellType();
    if(n_cell_type != 0 && n_cell_type != 4) continue;

	  
    uint32_t EID = essource_.emap_.detId2eid(tbRH_itr->id());
    HGCalTBElectronicsId eid(EID);
    int eSkiroc = (eid.iskiroc() - 1);

    //getting X and Y coordinates                                                                                                            
    CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots((tbRH_itr->id()).layer(), (tbRH_itr->id()).sensorIU(), (tbRH_itr->id()).sensorIV(), (tbRH_itr->id()).iu(), (tbRH_itr->id()).iv(), sensorsize);

    int TBlayer = tbRH_itr->id().layer();
    int layerID = getCellLayer(TBlayer);

    //fixme va flippato
    int cellID = getCellID(-CellCentreXY.second, CellCentreXY.first);

    std::pair<int, int> layerCell(layerID, cellID);
    //  std::cout << " start layer = " << tbRH_itr->id().layer() << " end layer = " << layerID << " TB id = " << eid.ichan() << " cellID = " << cellID << std::endl;
    
    unsigned int rawiddi = getRawID(layerCell);
    //    std::cout << " rawiddi = " << rawiddi << std::endl;
    DetId rechitID(rawiddi);

    //    std::cout << " build recHit "<< std::endl;
    HGCRecHit recHit( rechitID, tbRH_itr->energy(), tbRH_itr->time()); 
    float rhEnergy = ((tbRH_itr->energy() - commonmode[TBlayer-1]) / ADCtoMIP.at(eSkiroc) ) * (weights2GeV * weights2MIP * Weights_L.at(TBlayer-1) + 1. * MIP2GeV_sim);
    recHit.setEnergy(rhEnergy);	  

    //  std::cout << " push in event "<< std::endl;
    rechits->push_back(recHit);
  }
  event.put(std::move(rechits), outputCollectionName);
  
}
// Should there be a destructor ??
DEFINE_FWK_MODULE(HGCalRecHitProducer);
