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



void fitResolution(){
  gROOT->Reset();
  gROOT->Macro("~/public/setStyle.C");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);


  

  int nLayers = 28;
  int config = 1;
  int iColors[3] = {kRed, kBlue, kGreen+2};
  int energyP[3] = {20, 100, 200};



  /*    
  int nLayers = 8;
  int config = 1;
  int iColors[5] = {kRed, kCyan, kBlue, kGreen+2, kBlack};
  int energyP[5] = {20, 70, 100, 200, 250};
  */


  /*        
  int nLayers = 8;
  int config = 2;
  int iColors[5] = {kRed, kCyan, kBlue, kGreen+2, kBlack};
  int energyP[5] = {20, 70, 100, 200, 250};
  */


  TGraphErrors* tg[5];
  TGraphErrors* tgR[5];
  TGraphErrors* tgS[5];
  for(int iT=0; iT<5; ++iT){
    if(nLayers == 28 && iT >= 3) continue;
    tg[iT] = new TGraphErrors();
    tg[iT]->SetName(Form("mean_E%d", energyP[iT]));
    //    tg[iT]->SetPoint(0, -1, -1);

    //  tg[iT]->SetLineColor(iColors[iT]);
    //    tg[iT]->SetMarkerColor(iColors[iT]);
    tg[iT]->SetMarkerStyle(20);

    tgR[iT] = new TGraphErrors();
    tgR[iT]->SetName(Form("resolution_E%d", energyP[iT]));
    //    tgR[iT]->SetPoint(0, -1, -1);

    //    tgR[iT]->SetLineColor(iColors[iT]);
    //    tgR[iT]->SetMarkerColor(iColors[iT]);
    tgR[iT]->SetMarkerStyle(20);

    tgS[iT] = new TGraphErrors();
    tgS[iT]->SetName(Form("resolutionS_E%d", energyP[iT]));
    //    tgS[iT]->SetPoint(0, -1, -1);

    //    tgS[iT]->SetLineColor(iColors[iT]);
    //    tgS[iT]->SetMarkerColor(iColors[iT]);
    tgS[iT]->SetMarkerStyle(20);
  }


  TFile* inF[5];
  if(nLayers == 8 && config == 1){
    for(int iP=0; iP<5; ++iP)
      inF[iP] = TFile::Open(Form("energyLayer_ROOT/outF_Ele%d.root", energyP[iP]));
  }
  else if(nLayers == 8 && config == 2){
    for(int iP=0; iP<5; ++iP)
      inF[iP] = TFile::Open(Form("energyLayer_ROOT/outF_endShower_Ele%d.root", energyP[iP]));
  }
  else{
    for(int iP=0; iP<3; ++iP)
      inF[iP] = TFile::Open(Form("energyLayer_ROOT/outF_28layers_Ele%d.root", energyP[iP]));
  }

   TH1F* recoEnergy[5];
   TH1F* recoEnergyRel[5];

  TF1* fitF[5];
  for(int iP=0; iP<5; ++iP){
    if(nLayers == 28 && iP >=3) continue;
    fitF[iP] = new TF1(Form("E%d",energyP[iP]), "gaus", 0, 500);
  }


  TF1* fitFRel[5];
  for(int iP=0; iP<5; ++iP){
    if(nLayers == 28 && iP >=3) continue;
    fitFRel[iP] = new TF1(Form("E%dR",energyP[iP]), "gaus", 0, 500);
  }


  for(int iC=0; iC<5; ++iC){
    if(nLayers == 28 && iC >=3) continue;
    recoEnergy[iC] = (TH1F*)(inF[iC]->Get("recoEnergy"));
    recoEnergyRel[iC] = (TH1F*)(inF[iC]->Get("recoEnergyRel"));
    
    recoEnergy[iC]->SetLineColor(iColors[iC]);
    recoEnergy[iC]->SetLineWidth(2);

    recoEnergyRel[iC]->SetLineColor(iColors[iC]);
    recoEnergyRel[iC]->SetLineWidth(2);

    fitF[iC]->SetParameters(recoEnergy[iC]->Integral(), recoEnergy[iC]->GetMean(), recoEnergy[iC]->GetRMS());
    fitF[iC]->SetLineColor(iColors[iC]);

    fitFRel[iC]->SetParameters(recoEnergyRel[iC]->Integral(), recoEnergyRel[iC]->GetMean(), recoEnergyRel[iC]->GetRMS());
    fitFRel[iC]->SetLineColor(iColors[iC]);
  }


  std::cout << " ci sono " << std::endl;

  TLegend *legTGM = new TLegend(0.80,0.30,0.95,0.55,NULL,"brNDC");
  legTGM->SetTextFont(42);
  legTGM->SetFillColor(kWhite);
  legTGM->SetLineColor(kWhite);
  legTGM->SetShadowColor(kWhite);
  for(int iP=0; iP<5; ++iP){
    if(nLayers == 28 && iP >=3) continue;
    legTGM->AddEntry(recoEnergy[iP], Form("%dGeV", energyP[iP]), "l");
  }


  TLegend *leg = new TLegend(0.80,0.70,0.98,0.95,NULL,"brNDC");
  leg->SetTextFont(42);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  leg->SetShadowColor(kWhite);
  for(int iP=0; iP<5; ++iP){
    if(nLayers == 28 && iP >=3) continue;
    leg->AddEntry(recoEnergy[iP], Form("%dGeV", energyP[iP]), "l");
  }


  TLatex t1;
  t1.SetNDC();
  t1.SetTextSize(0.03);
  t1.SetTextFont(132);
  t1.SetTextColor(kRed);

  TLatex t2;
  t2.SetNDC();
  t2.SetTextSize(0.03);
  t2.SetTextFont(132);
  t2.SetTextColor(kBlue);


  TLatex t3;
  t3.SetNDC();
  t3.SetTextSize(0.03);
  t3.SetTextFont(132);
  t3.SetTextColor(kGreen+2);

  for(int iP=0; iP<5; ++iP){
    if(nLayers == 28 && iP >=3) continue;
    recoEnergy[iP]->Fit(Form("E%d", energyP[iP]), "R");
  }

  //    float energyValues[5] = {10., 30., 60., 100., 250.};

  //  std::string folder = "plots"+std::string(Form("_E%d",ene));
  std::string folder = "plots_Energy";
  //     std::string folder = "plots_E250";
  
  TCanvas* chE = new TCanvas();
  chE->cd();
  recoEnergy[0]->GetXaxis()->SetTitle("#Sigma E_{i}");
  recoEnergy[0]->Draw();
  for(int iP=1; iP<5; ++iP){
    if(nLayers == 28 && iP >=3) continue;
    recoEnergy[iP]->Draw("same");
  }
   t1.DrawLatex(0.2,0.9,Form("m = %.2e +/- %.2e   #sigma = %.2e +/- %.2e", fitF[0]->GetParameter(1), fitF[0]->GetParError(1), 
   			     fitF[0]->GetParameter(2), fitF[0]->GetParError(2)));
   t2.DrawLatex(0.2,0.85,Form("m = %.2e +/- %.2e   #sigma = %.2e +/- %.2e", fitF[1]->GetParameter(1), fitF[1]->GetParError(1), 
			       fitF[1]->GetParameter(2), fitF[1]->GetParError(2)));
   t3.DrawLatex(0.2,0.8,Form("m = %.2e +/- %.2e   #sigma = %.2e +/- %.2e", fitF[2]->GetParameter(1), fitF[2]->GetParError(1), 
			       fitF[2]->GetParameter(2), fitF[2]->GetParError(2)));

   /*
   for(int iT=0; iT<3; ++iT){
     tg[0]->SetPoint(iT+1, energyP[iT], fitF[iT]->GetParameter(1) / (1.* energyP[iT]));
     tg[0]->SetPointError(iT+1, 0., fitF[iT]->GetParError(1)/ (1.* energyP[iT]));

     tgR[0]->SetPoint(iT+1, 1./energyP[iT], pow(fitF[iT]->GetParameter(2),2) / (1.* energyP[iT]) * 100.);
     tgR[0]->SetPointError(iT+1, 0., 100.*pow(fitF[iT]->GetParError(2),2.)/ (1.* energyP[iT]));

     tgS[0]->SetPoint(iT+1, energyP[iT], fitF[iT]->GetParameter(2) / (1.* energyP[iT]));
     tgS[0]->SetPointError(iT+1, 0., fitF[iT]->GetParError(2) / (1.* energyP[iT]));


     //std::cout << " >>> energyP[iT] = " << energyP[iT] << " fitF[iT]->GetParameter(1) / (1.* energyP[iT]) = " << fitF[iT]->GetParameter(1) / (1.* energyP[iT]) << std::endl;
     //     std::cout << " >>> energyP[iT] = " << energyP[iT] << " Y = " << pow(fitF[iT]->GetParameter(2),2) / (1.* energyP[iT]) * 100. << std::endl;

     tg[0]->SetPoint(tg[iT]->GetN(), 400, 5);
     tgR[0]->SetPoint(tg[iT]->GetN(), 400, 5);
     tgS[0]->SetPoint(tg[iT]->GetN(), 400, 5);
  }
   */

   leg->Draw("same");
  //ch_100->Print((folder+"/h_EoP_hitTrack_100.png").c_str(), "png");
   chE->Print(Form((folder+"/energyPlots_layers%d_config%d.png").c_str(), nLayers, config), "png");


   /////////////////////////////
   
   for(int iP=0; iP<5; ++iP)
     for(int iP=0; iP<5; ++iP){
       if(nLayers == 28 && iP >=3) continue;
       recoEnergyRel[iP]->Fit(Form("E%dR", energyP[iP]), "R");
     }
  
   TCanvas* chER = new TCanvas();
   chER->cd();
   recoEnergyRel[0]->GetXaxis()->SetTitle("#Sigma E_{i} / E_{beam}");
   recoEnergyRel[0]->GetXaxis()->SetRangeUser(0.7, 1.3);
   recoEnergyRel[0]->GetYaxis()->SetRangeUser(0., 100.e+03);
   if(nLayers == 8){
     recoEnergyRel[0]->GetXaxis()->SetRangeUser(0., 1.1);
     recoEnergyRel[0]->GetYaxis()->SetRangeUser(0., 3.e+03);
   }
   recoEnergyRel[0]->Draw();
   for(int iP=1; iP<5; ++iP){
     if(nLayers == 28 && iP >=3) continue;
     recoEnergyRel[iP]->Draw("same");
   }
   
   /*
   t1.DrawLatex(0.2,0.9,Form("m = %.2e +/- %.2e   #sigma = %.2e +/- %.2e", fitFRel[0]->GetParameter(1), fitFRel[0]->GetParError(1), 
   			     fitFRel[0]->GetParameter(2), fitFRel[0]->GetParError(2)));
   t2.DrawLatex(0.2,0.85,Form("m = %.2e +/- %.2e   #sigma = %.2e +/- %.2e", fitFRel[1]->GetParameter(1), fitFRel[1]->GetParError(1), 
			       fitFRel[1]->GetParameter(2), fitFRel[1]->GetParError(2)));
   t3.DrawLatex(0.2,0.8,Form("m = %.2e +/- %.2e   #sigma = %.2e +/- %.2e", fitFRel[2]->GetParameter(1), fitFRel[2]->GetParError(1), 
			       fitFRel[2]->GetParameter(2), fitFRel[2]->GetParError(2)));
   */


   for(int iT=0; iT<5; ++iT){
     if(nLayers == 28 && iT >=3) continue;

     tg[0]->SetPoint(iT+1, energyP[iT], fitFRel[iT]->GetParameter(1));
     tg[0]->SetPointError(iT+1, 0., fitFRel[iT]->GetParError(1));

     tgR[0]->SetPoint(iT+1, 1./energyP[iT], pow(fitFRel[iT]->GetParameter(2),2) / (1.* energyP[iT]) * 100.);
     tgR[0]->SetPointError(iT+1, 0., 100.*pow(fitFRel[iT]->GetParError(2),2.)/ (1.* energyP[iT]));

     tgS[0]->SetPoint(iT+1, energyP[iT], fitFRel[iT]->GetParameter(2));
     tgS[0]->SetPointError(iT+1, 0., fitFRel[iT]->GetParError(2));


     //std::cout << " >>> energyP[iT] = " << energyP[iT] << " fitF[iT]->GetParameter(1) / (1.* energyP[iT]) = " << fitF[iT]->GetParameter(1) / (1.* energyP[iT]) << std::endl;
     //     std::cout << " >>> energyP[iT] = " << energyP[iT] << " Y = " << pow(fitF[iT]->GetParameter(2),2) / (1.* energyP[iT]) * 100. << std::endl;

     tg[0]->SetPoint(tg[iT]->GetN(), 400, 5);
     tgR[0]->SetPoint(tg[iT]->GetN(), 400, 5);
     tgS[0]->SetPoint(tg[iT]->GetN(), 400, 5);
  }
   leg->Draw("same");
  //ch_100->Print((folder+"/h_EoP_hitTrack_100.png").c_str(), "png");
   chER->Print(Form((folder+"/energyRelPlots_layers%d_config%d.png").c_str(), nLayers, config), "png");





    /////////TGraph
    TCanvas* tgM = new TCanvas();
    tgM->cd();
    tg[0]->GetYaxis()->SetTitle("#Sigma E_{i} / E_{beam}");
    tg[0]->GetXaxis()->SetTitle("E_{beam} (GeV)");
    tg[0]->GetXaxis()->SetRangeUser(0, 300.);
    tg[0]->GetYaxis()->SetRangeUser(0.99, 1.01);
    if(nLayers == 8) tg[0]->GetYaxis()->SetRangeUser(0., 1.01);
    tg[0]->Draw("ap");
    tgM->Print(Form((folder+"/energy_MeanVsEnergy_layers%d_config%d.png").c_str(), nLayers, config), "png");


    /*
    /////////TGraph
    TCanvas* tgReso = new TCanvas();
    tgReso->cd();
    tgR[0]->GetYaxis()->SetTitle("#sigma(E) / E)^{2} (%)");
    tgR[0]->GetXaxis()->SetTitle("#gamma 1/E (GeV)");
    tgR[0]->GetXaxis()->SetRangeUser(0, 0.15);
    tgR[0]->GetYaxis()->SetRangeUser(0., 1.);
    tgR[0]->Draw("ap");
    tgR[1]->Draw("p, same");
    tgR[2]->Draw("p, same");
    legTGM->Draw("same");
    tgReso->Print((folder+"/energy_ResolutionVsEnergy.png").c_str(), "png");
    */

    /////////TGraph   simple
    TCanvas* tgResoS = new TCanvas();
    tgResoS->cd();
    tgS[0]->GetYaxis()->SetTitle("#sigma( #Sigma E / E_{beam})");
    tgS[0]->GetXaxis()->SetTitle("E_{beam} (GeV)");
    tgS[0]->GetXaxis()->SetRangeUser(0., 300.);
    tgS[0]->GetYaxis()->SetRangeUser(0., 0.1);
    tgS[0]->Draw("ap");
    tgResoS->Print(Form((folder+"/energy_SResolutionVsEnergy_layers%d_config%d.png").c_str(), nLayers, config), "png");

    //    return;
    ////////////////
    TFile* inF2[5];
    if(nLayers == 28 && config == 1){
      for(int iP=0; iP<3; ++iP)
	inF2[iP] = TFile::Open(Form("showers_ROOT/analyzed_28layers_Ele%d.root", energyP[iP]));
    }
    else if(nLayers == 8 && config == 2){
      for(int iP=0; iP<5; ++iP)
	inF2[iP] = TFile::Open(Form("showers_ROOT/analyzed_endShower_Ele%d.root", energyP[iP]));
    }
    else{
      for(int iP=0; iP<5; ++iP)
	inF2[iP] = TFile::Open(Form("showers_ROOT/analyzed_Ele%d.root", energyP[iP]));
    }


    std::cout << " >>> fatto " << std::endl;
    TGraphErrors* tEF[5];
    for(int iT=0; iT<5; ++iT){
      if(nLayers == 28 && iT >= 3) continue;
      tEF[iT] = (TGraphErrors*)inF2[iT]->Get("beamEfraction_layer");
      tEF[iT]->SetName(Form("beamEfractionvsL_E%d", energyP[iT]));
      
      tEF[iT]->SetMarkerStyle(20);
      tEF[iT]->SetMarkerColor(iColors[iT]);
      tEF[iT]->SetLineColor(iColors[iT]);
      tEF[iT]->SetLineWidth(2);
      tEF[iT]->SetLineStyle(iT);
    }
    


    TCanvas* tgBeamEF = new TCanvas();
    tgBeamEF->cd();
    tEF[0]->GetXaxis()->SetTitle("shower depth (X_{0})");
    tEF[0]->GetYaxis()->SetTitle("#Sigma_{layer} E_{layer} / E_{beam}");
    tEF[0]->GetXaxis()->SetRangeUser(0., 300.);
    tEF[0]->GetYaxis()->SetRangeUser(0., 1.1);
    tEF[0]->Draw("ap");
    for(int iT=1; iT<5; ++iT){
      if(nLayers == 28 && iT >= 3) continue;
      tEF[iT]->Draw("p, same");
    }
    legTGM->Draw("same");
    tgBeamEF->Print(Form((folder+"/beamEnergyFractionLayer%d_config%d.png").c_str(), nLayers, config), "png");
    

}
