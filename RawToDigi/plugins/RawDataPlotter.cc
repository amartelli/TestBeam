#include <iostream>
#include "TH1F.h"
#include "TH2F.h"
#include "TH2Poly.h"
#include "TCanvas.h"
#include <TStyle.h>
#include <fstream>
#include <sstream>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "HGCal/DataFormats/interface/HGCalTBSkiroc2CMSCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#include <iomanip>
#include <set>


// struct hgcal_channel{
//   hgcal_channel(): Id(0),
// 		   iB(0),
// 		   iS(0),
// 		   iC(0),
//                    iAdcHG(0.),
//                    iAdcLG(0.){;}
//   int Id;
//   int IB, iS, iC;
//   std::vector<std::pair<std::pair<int, int>, float> > iAdcHG;
//   std::vector<std::pair<std::pair<int, int>, float> > iAdcLG;
// };

class RawDataPlotter : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  explicit RawDataPlotter(const edm::ParameterSet&);
  ~RawDataPlotter();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
private:
  virtual void beginJob() override;
  void analyze(const edm::Event& , const edm::EventSetup&) override;
  virtual void endJob() override;
  void InitTH2Poly(TH2Poly& poly, int layerID, int sensorIU, int sensorIV);

  struct {
    HGCalElectronicsMap emap_;
  } essource_;
  int m_sensorsize;
  bool m_eventPlotter;
  std::string m_electronicMap;

  int m_evtID;
  std::map<int,TH1F*> m_h_adcHigh;
  std::map<int,TH1F*> m_h_adcLow;
  std::map<int,TH1F*> m_h_adcHighDiff;
  std::map<int,TH1F*> m_h_adcLowDiff;
  std::map<int, std::map<std::pair<int, int>, std::vector<float> > > m_AdcHG;
  std::map<int, std::map<std::pair<int, int>, std::vector<float> > > m_AdcLG;
  std::map<int,TH2F*> m_ToAvsChannel;
  std::map<int,TH2F*> m_h2_adcHGVsSCA;
  std::map<int,TH2F*> m_h2_adcLGVsSCA;
  std::map<int,TH2F*> m_h2_adcHGVsSCAShift;
  std::map<int,TH2F*> m_h2_adcLGVsSCAShift;
  std::map<int,TH2F*> m_h_pulseHigh;
  std::map<int,TH2F*> m_h_pulseLow;

  edm::EDGetTokenT<HGCalTBSkiroc2CMSCollection> m_HGCalTBSkiroc2CMSCollection;

  HGCalTBTopology IsCellValid;
  HGCalTBCellVertices TheCell;
  std::vector<std::pair<double, double>> CellXY;
  std::pair<double, double> CellCentreXY;
  std::set< std::pair<int,HGCalTBDetId> > setOfConnectedDetId;
};

RawDataPlotter::RawDataPlotter(const edm::ParameterSet& iConfig) :
  m_sensorsize(iConfig.getUntrackedParameter<int>("SensorSize",128)),
  m_eventPlotter(iConfig.getUntrackedParameter<bool>("EventPlotter",false)),
  m_electronicMap(iConfig.getUntrackedParameter<std::string>("ElectronicMap","HGCal/CondObjects/data/map_CERN_Hexaboard_OneLayers_May2017.txt"))
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  m_HGCalTBSkiroc2CMSCollection = consumes<HGCalTBSkiroc2CMSCollection>(iConfig.getParameter<edm::InputTag>("InputCollection"));

  m_evtID=0;
  

  std::cout << " init iboard tot = " << HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD << std::endl;

  std::ostringstream os( std::ostringstream::ate );
  TH2F* htmp2;
  TH2F* htmp1vsSCA;
  TH1F* htmp1;
  //  TH1F* htmpValue;
  for(size_t ib = 0; ib<HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD; ib++) {
    //    if(ib > 0) return;
    for( size_t iski=0; iski<HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA; iski++ ){
      os.str("");os<<"HexaBoard"<<ib<<"_Skiroc"<<iski;
      TFileDirectory dir = fs->mkdir( os.str().c_str() );
 
      os.str("");
      os << "ToAvsChannel_HexaBoard" << ib << "_Chip" << iski ;
      htmp2=dir.make<TH2F>(os.str().c_str(),os.str().c_str(),128, 0., 128, 4096, 0., 4096.);
      m_ToAvsChannel.insert( std::pair<int,TH2F*>(ib*1000+iski*100, htmp2) );
      std::cout << " inizio ib = " << ib << " ski = " << iski << " vero iski = " << iski << std::endl;

     for( size_t ichan=0; ichan<HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC; ichan++ ){
	for( size_t it=0; it<NUMBER_OF_SCA; it++ ){
	  os.str("");
	  os << "HighGain_HexaBoard" << ib << "_Chip" << iski << "_Channel" << ichan << "_SCA" << it ;
	  htmp1=dir.make<TH1F>(os.str().c_str(),os.str().c_str(),1000,0,4096);
	  m_h_adcHigh.insert( std::pair<int,TH1F*>(ib*100000+iski*10000+ichan*100+it, htmp1) );
	  os.str("");
	  os << "LowGain_HexaBoard" << ib << "_Chip" << iski << "_Channel" << ichan << "_SCA" << it ;
	  htmp1=dir.make<TH1F>(os.str().c_str(),os.str().c_str(),1000,0,4096);
	  m_h_adcLow.insert( std::pair<int,TH1F*>(ib*100000+iski*10000+ichan*100+it, htmp1) );
	  
	  for(int ij=0; ij<13; ++ij){
	    std::pair<int, int> locPair(it, ij);
	    std::vector<float> adcV;
	    adcV.clear();
	    m_AdcHG[ib*1000+iski*100+ichan].insert(std::pair<std::pair<int, int>, std::vector<float> > (locPair, adcV));
	    m_AdcLG[ib*1000+iski*100+ichan].insert(std::pair<std::pair<int, int>, std::vector<float> > (locPair, adcV));

	    //if(iski==0 && ichan == 0 && (ib == 2 || ib == 3)){
	    /*
	    if(iski==0 && ichan == 0 && ib == 0){
	      os.str("");
	      os << "DiffHG_HexaBoard" << ib << "_Chip" << iski << "_Channel" << ichan << "_SCA" << it << "_TS" << ij;
	      htmpValue=dir.make<TH1F>(os.str().c_str(),os.str().c_str(),1000,0,4096);
	      m_h_adcHighDiff.insert( std::pair<int,TH1F*>(ib*1000+it*100+ij, htmpValue) );

	      os.str("");
	      os << "DiffLG_HexaBoard" << ib << "_Chip" << iski << "_Channel" << ichan << "_SCA" << it << "_TS" << ij;
	      htmpValue=dir.make<TH1F>(os.str().c_str(),os.str().c_str(),1000,0,4096);
	      m_h_adcLowDiff.insert( std::pair<int,TH1F*>(ib*1000+it*100+ij, htmpValue) );
	    }
	    */
	    
	  }
	}
	os.str("");
	os << "PedHG_TSvsSCA_HexaBoard" << ib << "_Chip" << iski << "_Channel" << ichan;
	htmp1vsSCA = dir.make<TH2F>(os.str().c_str(),os.str().c_str(),NUMBER_OF_SCA,0, NUMBER_OF_SCA, NUMBER_OF_SCA,0, NUMBER_OF_SCA);
	m_h2_adcHGVsSCA[ib*1000+iski*100+ichan] = htmp1vsSCA;
	os.str("");
	os << "PedLG_TSvsSCA_HexaBoard" << ib << "_Chip" << iski << "_Channel" << ichan;
	htmp1vsSCA = dir.make<TH2F>(os.str().c_str(),os.str().c_str(),NUMBER_OF_SCA,0, NUMBER_OF_SCA, NUMBER_OF_SCA,0, NUMBER_OF_SCA);
	m_h2_adcLGVsSCA[ib*1000+iski*100+ichan] = htmp1vsSCA;

	/*
	if(ib == 3){
	os.str("");
	os << "PedHG_TSvsSCA_HexaBoard_Shift" << ib << "_Chip" << iski << "_Channel" << ichan;
	htmp1vsSCA = dir.make<TH2F>(os.str().c_str(),os.str().c_str(),NUMBER_OF_SCA,0, NUMBER_OF_SCA, NUMBER_OF_SCA,0, NUMBER_OF_SCA);
	m_h2_adcHGVsSCAShift[ib*1000+iski*100+ichan] = htmp1vsSCA;
	os.str("");
	os << "PedLG_TSvsSCA_HexaBoard_Shift" << ib << "_Chip" << iski << "_Channel" << ichan;
	htmp1vsSCA = dir.make<TH2F>(os.str().c_str(),os.str().c_str(),NUMBER_OF_SCA,0, NUMBER_OF_SCA, NUMBER_OF_SCA,0, NUMBER_OF_SCA);
	m_h2_adcLGVsSCAShift[ib*1000+iski*100+ichan] = htmp1vsSCA;
	}
	*/

	if( ichan%2==0 ){
	  os.str("");
	  os << "HighGainVsSCA_Hexa" << ib << "_Chip" << iski << "_Channel" << ichan;
	  htmp2=dir.make<TH2F>(os.str().c_str(),os.str().c_str(),NUMBER_OF_SCA,0, NUMBER_OF_SCA,1000,0,4096);
	  m_h_pulseHigh.insert( std::pair<int,TH2F*>(ib*1000+iski*100+ichan, htmp2) );
	  os.str("");
	  os << "LowGainVsSCA_Hexa" << ib << "_Chip" << iski << "_Channel" << ichan;
	  htmp2=dir.make<TH2F>(os.str().c_str(),os.str().c_str(),NUMBER_OF_SCA,0, NUMBER_OF_SCA,1000,0,4096);
	  m_h_pulseLow.insert( std::pair<int,TH2F*>(ib*1000+iski*100+ichan, htmp2) );
	}
     }
    }
  }

  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(m_electronicMap);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  };
  std::cout << iConfig.dump() << std::endl;
}


RawDataPlotter::~RawDataPlotter()
{

}

void RawDataPlotter::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  edm::Handle<HGCalTBSkiroc2CMSCollection> skirocs;
  event.getByToken(m_HGCalTBSkiroc2CMSCollection, skirocs);

  if( !skirocs->size() ) return;
  
  std::map<int,TH2Poly*>  polyMap;
  if( m_eventPlotter ){
    std::ostringstream os( std::ostringstream::ate );
    os << "Event" << event.id().event();
    TFileDirectory dir = fs->mkdir( os.str().c_str() );
    for(size_t ib = 0; ib<HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD; ib++) {
      for( size_t it=0; it<NUMBER_OF_SCA; it++ ){
	TH2Poly *h=dir.make<TH2Poly>();
	os.str("");
	os<<"HexaBoard"<<ib<<"_SCA"<<it;
	h->SetName(os.str().c_str());
	h->SetTitle(os.str().c_str());
	InitTH2Poly(*h, (int)ib, 0, 0);
	polyMap.insert( std::pair<int,TH2Poly*>(100*ib+it,h) );
      }
    }
  }

  //  std::cout << " tot skiroc = " << skirocs->size() << " per board = " << HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA << std::endl;

  
  for( size_t iski=0;iski<skirocs->size(); iski++ ){
    HGCalTBSkiroc2CMS skiroc=skirocs->at(iski);
    std::vector<int> rollpositions=skiroc.rollPositions();
    int iboard=iski/HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;
    //    if(iboard > 0) return;
    for( size_t ichan=0; ichan<HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC; ichan++ ){
      HGCalTBDetId detid=skiroc.detid( ichan );
      HGCalTBElectronicsId eid( essource_.emap_.detId2eid(detid.rawId()) );
      if( essource_.emap_.existsEId(eid) ){
	std::pair<int,HGCalTBDetId> p( iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan,skiroc.detid(ichan) );
	setOfConnectedDetId.insert(p);
      }
      for( size_t it=0; it<NUMBER_OF_SCA; it++ ){
	if( rollpositions[it]<9 ){ //rm on track samples+2 last time sample which show weird behaviour
	  m_h_adcHigh[iboard*100000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*10000+ichan*100+it]->Fill(skiroc.ADCHigh(ichan,it));
	  m_h_adcLow[iboard*100000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*10000+ichan*100+it]->Fill(skiroc.ADCLow(ichan,it));
	}
	//	if(skiroc.TOAFall(ichan) == 4){
	  std::pair<int, int> locPair(it, rollpositions[it]);
	  m_AdcHG[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan][locPair].push_back(skiroc.ADCHigh(ichan,it));
	  m_AdcLG[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan][locPair].push_back(skiroc.ADCLow(ichan,it));
	  //	}
	  if( ichan%2==0 ){
	    m_h_pulseHigh[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan]->Fill( it,skiroc.ADCHigh(ichan,it) ); 
	    m_h_pulseLow[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan]->Fill( it,skiroc.ADCLow(ichan,it) );
	  }
	
	if(m_eventPlotter&&essource_.emap_.existsEId(eid) ){
	  if(!IsCellValid.iu_iv_valid( detid.layer(),detid.sensorIU(), detid.sensorIV(), detid.iu(), detid.iv(), m_sensorsize ) )  continue;
	  CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots( detid.layer(), detid.sensorIU(), detid.sensorIV(), detid.iu(), detid.iv(), m_sensorsize );
	  double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + HGCAL_TB_GEOMETRY::DELTA) : (CellCentreXY.first - HGCAL_TB_GEOMETRY::DELTA) ;
	  double iuy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + HGCAL_TB_GEOMETRY::DELTA) : (CellCentreXY.second - HGCAL_TB_GEOMETRY::DELTA);
	  polyMap[ 100*iboard+it ]->Fill(iux/2 , iuy, skiroc.ADCHigh(ichan,it) );
	}
      }
      //std::cout << "ib = " << iboard << " ski = " << iski << " ichan = " << ichan << " toa = " << skiroc.TOAFall(ichan) << std::endl;
      m_ToAvsChannel[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100]->Fill(ichan, skiroc.TOAFall(ichan));

    }
  }
}

void RawDataPlotter::InitTH2Poly(TH2Poly& poly, int layerID, int sensorIU, int sensorIV)
{
  double HexX[HGCAL_TB_GEOMETRY::MAXVERTICES] = {0.};
  double HexY[HGCAL_TB_GEOMETRY::MAXVERTICES] = {0.};
  for(int iv = -7; iv < 8; iv++) {
    for(int iu = -7; iu < 8; iu++) {
      if(!IsCellValid.iu_iv_valid(layerID, sensorIU, sensorIV, iu, iv, m_sensorsize)) continue;
      CellXY = TheCell.GetCellCoordinatesForPlots(layerID, sensorIU, sensorIV, iu, iv, m_sensorsize);
      assert(CellXY.size() == 4 || CellXY.size() == 6);
      unsigned int iVertex = 0;
      for(std::vector<std::pair<double, double>>::const_iterator it = CellXY.begin(); it != CellXY.end(); it++) {
	HexX[iVertex] =  it->first;
	HexY[iVertex] =  it->second;
	++iVertex;
      }
      //Somehow cloning of the TH2Poly was not working. Need to look at it. Currently physically booked another one.
      poly.AddBin(CellXY.size(), HexX, HexY);
    }//loop over iu
  }//loop over iv
}

void RawDataPlotter::beginJob()
{
}

void RawDataPlotter::endJob()
{
  for(size_t iboard = 0; iboard<HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD; iboard++) {
    //    if(iboard > 0) return;
    for( size_t iski=0; iski<HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA; iski++ ){
      for(size_t ichan = 0; ichan<HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC; ichan++ ){
	for( size_t it=0; it<NUMBER_OF_SCA; it++ ){
	  // int sizeLoc = m_AdcHG[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan].size();
	  // std::cout << " sizeLoc = " << sizeLoc << std::endl;
	  for(int ij=0; ij<13; ++ij){
	    std::pair<int, int> locPair(it, ij);
	    std::vector<float> HGv(m_AdcHG[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan][locPair]);

	    if(HGv.size() == 0){
	      m_h2_adcHGVsSCA[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan]->Fill(it, ij, -1);
	      // if(iboard == 3 && ij < 4) m_h2_adcHGVsSCAShift[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan]->Fill(it, ij+9, -1);
	      // if(iboard == 3 && ij >= 4) m_h2_adcHGVsSCAShift[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan]->Fill(it, ij-4, -1);
	    }
	    else{
	      std::sort(HGv.begin(), HGv.end());
	      m_h2_adcHGVsSCA[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan]->Fill(it, ij, HGv.at(1.*HGv.size()/2.));
	      // if(iboard == 3 && ij < 4) m_h2_adcHGVsSCAShift[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan]->Fill(it, ij+9, HGv.at(1.*HGv.size()/2.));
	      // if(iboard == 3 && ij >= 4) m_h2_adcHGVsSCAShift[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan]->Fill(it, ij-4, HGv.at(1.*HGv.size()/2.));
	    }

	    std::vector<float> LGv(m_AdcLG[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan][locPair]);
	    if(LGv.size() == 0){
	      m_h2_adcLGVsSCA[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan]->Fill(it, ij, -1);
	      // if(iboard == 3 && ij < 4) m_h2_adcLGVsSCAShift[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan]->Fill(it, ij+9, -1);
	      // if(iboard == 3 && ij >= 4) m_h2_adcLGVsSCAShift[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan]->Fill(it, ij-4, -1);
	    }
	    else{	    
	      std::sort(LGv.begin(), LGv.end());
	      m_h2_adcLGVsSCA[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan]->Fill(it, ij, LGv.at(1.*LGv.size()/2.));
	      // if(iboard == 3 && ij < 4) m_h2_adcLGVsSCAShift[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan]->Fill(it, ij+9, LGv.at(1.*LGv.size()/2.));
	      // if(iboard == 3 && ij >= 4) m_h2_adcLGVsSCAShift[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan]->Fill(it, ij-4, LGv.at(1.*LGv.size()/2.));
	    }

	    
	    //  if(ichan == 0 && iski == 0 && (iboard == 2 || iboard == 3) ){
	    /*
	    if(ichan == 0 && iski == 0 && iboard == 0){
	      //std::cout << " size HG = " << HGv.size() << std::endl;
	      for(unsigned int jk=0; jk<HGv.size(); ++jk) m_h_adcHighDiff[iboard*1000+it*100+ij]->Fill(HGv.at(jk));
	      for(unsigned int jk=0; jk<LGv.size(); ++jk) m_h_adcLowDiff[iboard*1000+it*100+ij]->Fill(LGv.at(jk));
	      //	      std::cout << " post size HG = " << HGv.size() << std::endl;
	      gStyle->SetOptStat(1);
	      TCanvas* cHisto = new TCanvas();
	      cHisto->cd();
	      m_h_adcHighDiff[iboard*1000+it*100+ij]->GetXaxis()->SetTitle("ADC HG");
	      m_h_adcHighDiff[iboard*1000+it*100+ij]->Draw();
	      cHisto->Print(Form("folderRUN/histos/%s.png", m_h_adcHighDiff[iboard*1000+it*100+ij]->GetName() ), ".png");
	      
	      TCanvas* cListo = new TCanvas();
	      cListo->cd();
	      m_h_adcLowDiff[iboard*1000+it*100+ij]->GetXaxis()->SetTitle("ADC LG");
	      m_h_adcLowDiff[iboard*1000+it*100+ij]->Draw();
	      cListo->Print(Form("folderRUN/histos/%s.png", m_h_adcLowDiff[iboard*1000+it*100+ij]->GetName() ), ".png");
	    }
	    */
	    
	  }
	}
	
	gStyle->SetOptStat(0);
	
	TCanvas* cH = new TCanvas();
	cH->cd();
	m_h2_adcHGVsSCA[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan]->GetXaxis()->SetTitle("SCA");
	m_h2_adcHGVsSCA[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan]->GetYaxis()->SetTitle("time stamp");
	m_h2_adcHGVsSCA[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan]->Draw("coltext");
	cH->Print(Form("folderRUN/%s.png", m_h2_adcHGVsSCA[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan]->GetName() ), ".png");

	TCanvas* cL = new TCanvas();
	cL->cd();
	m_h2_adcLGVsSCA[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan]->GetXaxis()->SetTitle("SCA");
	m_h2_adcLGVsSCA[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan]->GetYaxis()->SetTitle("time stamp");
	m_h2_adcLGVsSCA[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan]->Draw("coltext");
	cL->Print(Form("folderRUN/%s.png", m_h2_adcLGVsSCA[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan]->GetName() ), ".png");
	/*
	if(iboard == 3){
	  TCanvas* cHS = new TCanvas();
	  cHS->cd();
	  m_h2_adcHGVsSCAShift[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan]->GetXaxis()->SetTitle("SCA");
	  m_h2_adcHGVsSCAShift[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan]->GetYaxis()->SetTitle("time stamp - 4");
	  m_h2_adcHGVsSCAShift[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan]->Draw("coltext");
	  cHS->Print(Form("folderRUN/Shift/%s.png", m_h2_adcHGVsSCAShift[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan]->GetName() ), ".png");
	  
	  TCanvas* cLS = new TCanvas();
	  cLS->cd();
	  m_h2_adcLGVsSCAShift[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan]->GetXaxis()->SetTitle("SCA");
	  m_h2_adcLGVsSCAShift[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan]->GetYaxis()->SetTitle("time stamp - 4");
	  m_h2_adcLGVsSCAShift[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan]->Draw("coltext");
	  cLS->Print(Form("folderRUN/Shift/%s.png", m_h2_adcLGVsSCA[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan]->GetName() ), ".png");
	}
	*/
	
      }
      std::cout << " stampo ib = " << iboard << " ski = " << iski << " vero iski = " << iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA << std::endl;
      TCanvas* cToA = new TCanvas();
      cToA->cd();
      gPad->SetLogy();
      m_ToAvsChannel[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100]->GetXaxis()->SetTitle("channel");
      m_ToAvsChannel[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100]->GetYaxis()->SetTitle("ToA");
      m_ToAvsChannel[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100]->Draw("colz");
      cToA->Print(Form("folderRUN/TOA/%s.png", m_ToAvsChannel[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100]->GetName() ), ".png");
      //      cToA->Print(Form("folderRUN/TOA/%s.root", m_ToAvsChannel[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100]->GetName() ), ".root");
    }
  }
}

void RawDataPlotter::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(RawDataPlotter);
