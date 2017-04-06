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
  _layers_config(cfg.getParameter<int>("layers_config")){ 
  //  _treeName(cfg.getParameter<std::string>("treeName")){

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
      //      CMSSW_cellLayer = cfg.getParameter<std::vector<int> >("CMSSW_cellLayer_1");
      Weights_L = cfg.getParameter<std::vector<double> >("LayerWeight_8L_conf1");
    }
    if(_layers_config == 2) {
      //      CMSSW_cellLayer = cfg.getParameter<std::vector<int> >("CMSSW_cellLayer_2");
      Weights_L = cfg.getParameter<std::vector<double> >("LayerWeight_8L_conf2");
    }
  }
  //  std::cout << " HGCalRecHitProducer::HGCalRecHitProducer get File " << std::endl;

  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(_mapFile);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  };

    TBcmsswGeometryMap = new RecHitGeometryMapping(cfg);  

    /*
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  // h_layer = fs->make<TH1F>("h_layer", "", 10, 0., 10.);
  // h_IU = fs->make<TH1F>("h_IU", "", 100, 0., 100.);
  // h_IV = fs->make<TH1F>("h_IV", "", 100, 0., 100.);
  // h_iu = fs->make<TH1F>("h_iu", "", 100, 0., 100.);
  // h_iv = fs->make<TH1F>("h_iv", "", 100, 0., 100.);
  // h_cellType = fs->make<TH1F>("h_cellType", "", 10, 0., 10.);


  TBtree = fs->make<TTree>("TBtree", "TBtree");
  TBtree->Branch("layerTB",&layerTB, "layerTB/I");
  TBtree->Branch("localXTB",&localXTB, "localXTB/F");
  TBtree->Branch("localYTB",&localYTB, "localYTB/F");
  TBtree->Branch("detIDTB",&detIDTB, "detIDTB/i");
    */
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
  // float maxEnergy[8] = {-999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999., -999.};
  // float XmaxRH[8] = {-999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999., -999.};
  // float YmaxRH[8] = {-999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999., -999.}; 

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

    // h_layer->Fill((tbRH_itr->id()).layer());
    // h_IU->Fill((tbRH_itr->id()).sensorIU());    
    // h_IV->Fill((tbRH_itr->id()).sensorIV());
    // h_iu->Fill((tbRH_itr->id()).iu());    
    // h_iv->Fill((tbRH_itr->id()).iv());
    // h_cellType->Fill(n_cell_type);


    if(n_cell_type != 0 && n_cell_type != 4) continue;

	  
    uint32_t EID = essource_.emap_.detId2eid(tbRH_itr->id());
    HGCalTBElectronicsId eid(EID);
    int eSkiroc = (eid.iskiroc() - 1);

    // //getting X and Y coordinates                                                                                                            
    // CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots((tbRH_itr->id()).layer(), (tbRH_itr->id()).sensorIU(), (tbRH_itr->id()).sensorIV(), (tbRH_itr->id()).iu(), (tbRH_itr->id()).iv(), sensorsize);    

    int TBlayer = tbRH_itr->id().layer();

    unsigned int rawiddi = TBcmsswGeometryMap->getCMSSWID_fromTBID(tbRH_itr->id().rawId());
    DetId rechitID(rawiddi);

    //FIXME
    //DetId rechitID(tbRH_itr->id().rawId());

    //    std::cout << " build recHit "<< std::endl;
    HGCRecHit recHit( rechitID, tbRH_itr->energy(), tbRH_itr->time()); 
    float rhEnergy = ((tbRH_itr->energy() - commonmode[TBlayer-1]) / ADCtoMIP.at(eSkiroc) ) * (weights2GeV * weights2MIP * Weights_L.at(TBlayer-1) + 1. * MIP2GeV_sim);
    recHit.setEnergy(rhEnergy);	  

    //  std::cout << " push in event "<< std::endl;
    rechits->push_back(recHit);

  }

  // for(int ij=0; ij<8; ++ij)
  //   std::cout << "layer = " << ij+1 << " energy = " << maxEnergy[ij] << " XmaxRH = " << XmaxRH[ij] << " YmaxRH = " << YmaxRH[ij] << std::endl;


  event.put(std::move(rechits), outputCollectionName);

  /*  
  for(int iLay = 1; iLay<9; ++iLay){
    for(int jiu=-7; jiu<8; ++jiu){
      for(int jiv=-7; jiv<8; ++jiv){
	for(int kT=0; kT<6; ++kT){

	  HGCalTBDetId localDetID(iLay, 0, 0, jiu, jiv, kT);
	  CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots(iLay, 0, 0, jiu, jiv, sensorsize);       
  
	  layerTB = iLay;
	  localXTB = CellCentreXY.second;
	  localYTB = -CellCentreXY.first;
	  detIDTB = localDetID.rawId();
	  if(localXTB > -7 && localXTB < 7 && localYTB > -7 && localYTB < 7)  TBtree->Fill();
	}
      }
    }
  }
  */
}
// Should there be a destructor ??
DEFINE_FWK_MODULE(HGCalRecHitProducer);
