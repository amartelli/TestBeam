#include "TLegend.h"
#include "TLatex.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include "TTree.h"
#include "TChain.h"
#include <vector>
#include <fstream>
#include <string>
#include "TROOT.h"
#include "TSystem.h"

#include "/afs/cern.ch/work/a/amartell/ES_P5/CMSSW_7_4_6_patch6/src/MyAnalyzer/ESAnalyzer/macroAnalysis/langausfit.h"

void doPlots(int nLayer, int energyEle, int config){
  gROOT->Reset();
  gROOT->Macro("~/public/setStyle.C");
  gStyle->SetOptTitle(0);


  float dEdX_weights[28];
  float X0val[28];

  if(nLayer == 28){
    //    float dEdX_weights_28[28] = {7.28, 7.28, 7.28, 7.28, 7.28, 7.28, 7.28, 7.28, 7.28, 8.164, 9.347, 9.347, 9.347, 9.347, 9.347, 9.347, 9.347, 9.347, 9.347, 10.894, 12.866, 12.866, 12.866, 12.866, 12.866, 12.866, 12.866, 8.167};
    
    float dEdX_weights_28[28] = {7.203, 7.203, 7.203, 7.203, 7.203, 7.203, 7.203, 7.203, 7.203, 8.087, 9.27, 9.27, 9.27, 9.27, 9.27, 9.27, 9.27, 9.27, 9.27, 10.817, 12.789, 12.789, 12.789, 12.789, 12.789, 12.789, 12.789, 8.148};
    //float dEdX_weights_28[28] = {7.164/2., 7.164, 7.164, 7.164, 7.164, 7.164, 7.164, 7.164, 7.164, 8.048, 9.231, 9.231, 9.231, 9.231, 9.231, 9.231, 9.231, 9.231, 9.231, 10.778, 12.75, 12.75, 12.75, 12.75, 12.75, 12.75, 12.75, 8.109};
    float X0val_28[28] = {0.574, 0.699, 0.574, 0.699, 0.574, 0.699, 0.574, 0.699, 0.574, 0.699, 0.802, 0.977, 0.802, 0.977, 0.802, 0.977, 0.802, 0.977, 0.802, 0.977, 1.202, 1.439, 1.202, 1.439, 1.202, 1.439, 1.202, 1.439};

    for(int j=0; j<28; ++j){
      dEdX_weights[j] = dEdX_weights_28[j];
      X0val[j] = X0val_28[j];
    }
  }
  else{
    if(config == 1){
      //    float dEdX_weights_8[28] = {33.151, 13.261, 14.247, 9.885, 9.844, 9.844, 16.416, 28.296, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
      float dEdX_weights_8[28] = {33.074, 13.184, 14.17, 9.788, 9.766, 9.766, 16.339, 14.129, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
      float X0val_8[28] = {6.268, 1.131, 1.131, 1.362, 0.574, 1.301, 0.574, 2.42, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
      
      for(int j=0; j<28; ++j){
	dEdX_weights[j] = dEdX_weights_8[j];
	X0val[j] = X0val_8[j];
      }
    }
    else if(config == 2){
     float dEdX_weights_8[28] = {35.866, 30.864, 28.803, 23.095, 20.657, 19.804, 36.322, 27.451, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
     float X0val_8[28] = {5.048, 3.412, 3.412, 2.866, 2.512, 1.625, 2.368, 6.021, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
      
      for(int j=0; j<28; ++j){
	dEdX_weights[j] = dEdX_weights_8[j];
	X0val[j] = X0val_8[j];
      }
    }
  }

  //  float PDGMIP2GeV_sim = 5.28e-05;
  //  float MIP2GeV_sim = 5.591e-05;

  float MIP2GeV_sim = 5.28e-05;             // MIP in MC
  float weights2GeV = 1.e-03;
  //  float weights2MIP = 52.8 / 77.52;         // weights with MIP from PDG
  float weights2MIP = 52.8 / 63.6;         // weights with MIP from PDG
  //  float MIP2GeV_sim = 64.19e-06;
  //  float MIP2GeV_sim = 53.16e-06;
  //  float MIP2GeV_sim = 53.16e-06;


  //  weights2GeV = weights2GeV / 0.828;
  float G4Escale = 1.; // 0.83;   // >>> fix scale in data mc for Si
  //float G4Escale = 1.214;



  TChain* t = new TChain("HGCalTBAnalyzer/HGCTB");
  
  if(nLayer == 8 && config == 1){
    std::string nameFolder = "160926_133629";
    if(energyEle == 200) nameFolder = "160926_133629";
    if(energyEle == 20) nameFolder = "160926_133600";
    if(energyEle == 250) nameFolder = "160926_133635";
    if(energyEle == 70) nameFolder = "160926_133612";
    if(energyEle == 100) nameFolder = "160926_133618";
    
    for(int iFo=0; iFo<2; ++iFo){
      for(int iF=0; iF<1000; ++iF){
	t->Add(Form(("root://eoscms.cern.ch///store/group/upgrade/HGCAL/simulation/8moduleIv1/mc/CRAB_PrivateMC/crab_Ele%dGeV/"+nameFolder+"/000%d/TBGenSim_%d.root").c_str(), energyEle, iFo, iF));
      }
    }
  }
  else if(nLayer == 8 && config == 2){
    std::string nameFolder = "160927_101523";
    if(energyEle == 200) nameFolder = "160927_101553";
    if(energyEle == 20) nameFolder = "160927_101523";
    if(energyEle == 250) nameFolder = "160927_101559";
    if(energyEle == 70) nameFolder = "160927_101535";
    if(energyEle == 100) nameFolder = "160927_101541";

    for(int iFo=0; iFo<2; ++iFo){
      for(int iF=0; iF<1000; ++iF){
	t->Add(Form(("root://eoscms.cern.ch///store/group/upgrade/HGCAL/simulation/8moduleIIv1/mc/CRAB_PrivateMC/crab_Ele%dGeV/"+nameFolder+"/000%d/TBGenSim_%d.root").c_str(), energyEle, iFo, iF));
      }
    }    
  }
  else{
    std::string nameFolder = "161004_152035";
    if(energyEle == 100) nameFolder = "161004_152035";
    if(energyEle == 200) nameFolder = "161004_152044";
    if(energyEle == 20) nameFolder = "161004_152027";

    for(int iFo=0; iFo<2; ++iFo){
      for(int iF=0; iF<1000; ++iF){
	t->Add(Form(("root://eoscms.cern.ch///store/group/upgrade/HGCAL/simulation/28modulev1/mc/CRAB_PrivateMC/crab_Ele%dGeV/"+nameFolder+"/000%d/TBGenSim_%d.root").c_str(), energyEle, iFo, iF));
      }
    }
  }
  /*
  if(1 == 2){
    if(energyEle == 100){
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele100GeV/160831_132138/0000/TBGenSim_391.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele100GeV/160831_132138/0000/TBGenSim_448.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele100GeV/160831_132138/0000/TBGenSim_451.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele100GeV/160831_132138/0000/TBGenSim_459.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele100GeV/160831_132138/0000/TBGenSim_462.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele100GeV/160831_141056/0000/TBGenSim_1.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele100GeV/160831_141056/0000/TBGenSim_2.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele100GeV/160831_141056/0000/TBGenSim_3.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele100GeV/160831_141056/0000/TBGenSim_4.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele100GeV/160831_141056/0000/TBGenSim_5.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele100GeV/160831_141056/0000/TBGenSim_6.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele100GeV/160831_141056/0000/TBGenSim_7.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele100GeV/160831_141056/0000/TBGenSim_10.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele100GeV/160831_142922/0000/TBGenSim_1.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele100GeV/160831_142922/0000/TBGenSim_2.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele100GeV/160831_142922/0000/TBGenSim_3.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele100GeV/160831_142922/0000/TBGenSim_4.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele100GeV/160831_142922/0000/TBGenSim_5.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele100GeV/160831_142922/0000/TBGenSim_6.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele100GeV/160831_142922/0000/TBGenSim_7.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele100GeV/160831_142922/0000/TBGenSim_8.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele100GeV/160831_142922/0000/TBGenSim_9.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele100GeV/160831_142922/0000/TBGenSim_10.root");
    }
    if(energyEle == 20){
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele20GeV/160903_133956/0000/TBGenSim_1.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele20GeV/160903_133956/0000/TBGenSim_2.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele20GeV/160903_133956/0000/TBGenSim_3.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele20GeV/160903_133956/0000/TBGenSim_4.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele20GeV/160903_133956/0000/TBGenSim_5.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele20GeV/160903_133956/0000/TBGenSim_6.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele20GeV/160903_133956/0000/TBGenSim_7.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele20GeV/160903_133956/0000/TBGenSim_8.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele20GeV/160903_133956/0000/TBGenSim_9.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele20GeV/160903_133956/0000/TBGenSim_10.root");
    }
    if(energyEle == 250){
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele250GeV/160903_134005/0000/TBGenSim_1.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele250GeV/160903_134005/0000/TBGenSim_2.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele250GeV/160903_134005/0000/TBGenSim_3.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele250GeV/160903_134005/0000/TBGenSim_4.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele250GeV/160903_134005/0000/TBGenSim_5.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele250GeV/160903_134005/0000/TBGenSim_6.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele250GeV/160903_134005/0000/TBGenSim_7.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele250GeV/160903_134005/0000/TBGenSim_8.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele250GeV/160903_134005/0000/TBGenSim_9.root");
      t->Add("root://eoscms.cern.ch//store/user/shilpi/HGCALTB/28module_31Aug/mc/CRAB_PrivateMC/crab_Ele250GeV/160903_134005/0000/TBGenSim_10.root");
    }
  }
  */

  double totEntries = t->GetEntries();
  std::cout << " tree read => entries = " << totEntries << std::endl;

  TH1F* enLayer[28];
  TH1F* enLayerCorr[28];
  for(int i=0; i<28; ++i){
    if(nLayer != 28 && i >= 8) continue;
    enLayer[i] = new TH1F(Form("enLayer%d", i+1), "", 5000, 0., 0.5);
    enLayerCorr[i] = new TH1F(Form("enLayerCorr%d", i+1), "", 5000, 0., 10.);
  }

  TH1F* recoEnergy = new TH1F("recoEnergy", "", 5000, 0., 500.);
  TH1F* recoEnergyRel = new TH1F("recoEnergyRel", "", 1000, 0., 2.);

  double xBeam, yBeam;
  std::vector<float>* simHitLayEn2E = 0;
  t->SetBranchAddress("simHitLayEn2E", &simHitLayEn2E);
  t->SetBranchAddress("xBeam", &xBeam);
  t->SetBranchAddress("yBeam", &yBeam);


  //looping over entries
  for(int iE=0; iE<totEntries; ++iE){

    t->GetEntry(iE);
    //    std::cout << " simHitLayEn2E->size() = " << simHitLayEn2E->size() << std::endl;
    if(abs(xBeam) >= 2 || abs(yBeam) >= 2) continue;
    float totEnergy = 0.;
    for(int iS=0; iS<simHitLayEn2E->size(); ++iS){
      enLayer[iS]->Fill(simHitLayEn2E->at(iS));
      
      //      std::cout << " E = " simHitLayEn2E->at(iS) << << std::endl;
      float localE = simHitLayEn2E->at(iS) / G4Escale * ( 1./MIP2GeV_sim * weights2GeV * weights2MIP* dEdX_weights[iS] + 1.);

      enLayerCorr[iS]->Fill(localE);
      totEnergy += localE;
    }
    recoEnergy->Fill(totEnergy);
    recoEnergyRel->Fill(totEnergy/(1.*energyEle));
  }


  TFile* outF;
  if(nLayer == 28) outF = new TFile(Form("energyLayer_ROOT/outF_28layers_Ele%d.root", energyEle), "recreate");
  else if(config == 1) outF = new TFile(Form("energyLayer_ROOT/outF_Ele%d.root", energyEle), "recreate");
  else if(config == 2) outF = new TFile(Form("energyLayer_ROOT/outF_endShower_Ele%d.root", energyEle), "recreate");
  outF->cd();
  for(int i=0; i<28; ++i){
    if(nLayer != 28 && i >= 8) continue;
    enLayer[i]->Write();
    enLayerCorr[i]->Write();
  }
  recoEnergy->Write();
  recoEnergyRel->Write();
  outF->Close();


}
