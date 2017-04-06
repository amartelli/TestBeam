#include "HGCal/Reco/interface/RecHitGeometryMapping.h"
//#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
//#include "DataFormats/plugins/HGCGeometryReader.cc"
#include <iostream>

#include "TTree.h"
#include "TFile.h"
#include "TH2F.h"

using namespace std;

RecHitGeometryMapping::RecHitGeometryMapping(const edm::ParameterSet& cfg): 
  _layers_config(cfg.getParameter<int>("layers_config")), 
  _treeName(cfg.getParameter<std::string>("treeName")),
  _treeNameTB(cfg.getParameter<std::string>("treeNameTB")){
    
  //fixme check order 
  if(_layers_config == 1){
    CMSSW_cellLayer = cfg.getParameter<std::vector<int> >("CMSSW_cellLayer_1");
    //    Weights_L = cfg.getParameter<std::vector<double> >("LayerWeight_8L_conf1");
  }
  if(_layers_config == 2) {
    CMSSW_cellLayer = cfg.getParameter<std::vector<int> >("CMSSW_cellLayer_2");
    //    Weights_L = cfg.getParameter<std::vector<double> >("LayerWeight_8L_conf2");
  }

  //  std::cout << " RecHitGeometryMapping::RecHitGeometryMapping get File " << std::endl;

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

  /////// TB info
  TFile* inFTB = TFile::Open(_treeNameTB.c_str());
  TTree* treeTB = (TTree*)inFTB->Get("HGCalRecHit/TBtree");
  //  std::cout << " entries = " << tree->GetEntries() << std::endl;                                                                                                                                                                            
  Int_t   layerTB;
  Float_t localXTB, localYTB;
  UInt_t detIDTB;

  treeTB->SetBranchAddress("layerTB",&layerTB);
  treeTB->SetBranchAddress("detIDTB", &detIDTB);
  treeTB->SetBranchAddress("localXTB",&localXTB);
  treeTB->SetBranchAddress("localYTB",&localYTB);


  //  waferMap = new TH2F("waferMap", "",  200, -10, 10, 200, -10, 10);
  for(int ij=0; ij<tree->GetEntries(); ++ij){
    tree->GetEntry(ij);
    if(wafer == 426 && z > 0 && subdet == 3 && cell < 132){    
      if(layer > 0 && layer < 30){

	for(int jk=0; jk<treeTB->GetEntries(); ++jk){
	  treeTB->GetEntry(jk);

	  if(locX == localXTB && locY == localYTB && layerTB == getTBLayer(layer)){
	    CMSSWIDk_TBIDv[detID] = detIDTB;
	    break;
	  }
	}

      }
    }
  }


}

  
int RecHitGeometryMapping::getTBLayer(int layer){
  int TBlayer = 0;
  for(unsigned int iVal=0; iVal<CMSSW_cellLayer.size(); ++iVal){
    if(CMSSW_cellLayer.at(iVal) == layer) {
      TBlayer = iVal+1;
      break;
    }
  }
  
  return TBlayer; 
}


int RecHitGeometryMapping::getCMSSWLayer(int TBlayer){
  //std::cout << " TBlayer = " << TBlayer << "corresponds to = " << CMSSW_cellLayer.at(TBlayer-1) << std::endl;                                                                                                                         
  return CMSSW_cellLayer.at(TBlayer-1);
}


unsigned int RecHitGeometryMapping::getTBID_fromCMSSWID(unsigned int rawiddi){
  return CMSSWIDk_TBIDv[rawiddi];
  
}

unsigned int RecHitGeometryMapping::getCMSSWID_fromTBID(unsigned int rawiddiTB){
  unsigned int CMSSWID = 0;
  for(std::map<unsigned int, unsigned int>::const_iterator itr = CMSSWIDk_TBIDv.begin(); itr != CMSSWIDk_TBIDv.end(); ++itr){
    if(itr->second == rawiddiTB){
      CMSSWID = itr->first; 
      break;
    }
  }
  return CMSSWID;
}


// Should there be a destructor ??
//DEFINE_FWK_MODULE(RecHitGeometryMapping);
