#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>


#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TRandom3.h"


#include "interface/ZGSample.h"
#include "interface/ZGConfig.h"
#include "interface/ZGCommonTools.h"

// muon rochester corrections:
#include "../interface/rochcor2015.h"
#include "../interface/muresolution_run2.h"

// photon smearing/scales for 76X
#include "../interface/EnergyScaleCorrection_class.hh"


#define zg_cxx
#include "interface/zg.h"


#define DATABLINDING false

UInt_t DEBUG_EVENT = 1276086;

bool doSyst = false;


void addTreeToFile( TFile* file, ZGSample sample, const ZGConfig& cfg );
void addTreeToFile( TFile* file, const std::string& treeName, std::vector<ZGSample> samples, const ZGConfig& cfg, int idMin=-1, int idMax=-1 );
void smearEmEnergy     ( TLorentzVector& p );
void applyEmEnergyScale( TLorentzVector& p );
float getPUweight( TH1D* h1_pu, float nTrueInt );

//TLorentzVector selectPhoton( const ZGConfig& cfg, const ZGTree& myTree, int index, const TLorentzVector& lept0, const TLorentzVector& lept1 );
TLorentzVector selectPhoton( const ZGConfig& cfg, const ZGTree& myTree, int index, const TLorentzVector& lept0, const TLorentzVector& lept1, EnergyScaleCorrection_class egcor );


TRandom3 myRandom_(13);


int main( int argc, char* argv[] ) {


  std::cout << std::endl << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|               Running runZGAnalysis                |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "|                                                    |" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;




  if( argc<2 ) {
    std::cout << "USAGE: ./runZGAnalysis [configFileName] [data/MC]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  std::string configFileName(argv[1]);
  ZGConfig cfg(configFileName);



  bool onlyData = false;
  bool onlyMC   = false;
  bool onlySignal = false;
  bool noSignals = false;
  if( argc > 2 ) {

    std::string dataMC(argv[2]);
    if( dataMC=="data" ) onlyData = true;
    else if( dataMC=="MC" || dataMC=="mc" ) onlyMC = true;
    else if( dataMC=="mcbg" || dataMC=="mcBG" || dataMC=="MCBG" || dataMC=="mc_bg" || dataMC=="MC_BG" ) {
      onlyMC = true;
      noSignals = true;
    } else if( dataMC=="signal" ) onlySignal = true;
    else {
      std::cout << "-> You passed a second argument that isn't 'data', nor 'MC', nor 'signal', so I don't know what to do about it." << std::endl;
    }

  } else {

    std::cout << "-> Will run on both data and MC." << std::endl;

  }


  if( onlyMC ) {
    std::cout << "-> Will run only on MC." << std::endl;
    if( noSignals ) {
      std::cout << "-> Will skip signal." << std::endl;
    }
  } 

  if( onlyData ) {
    std::cout << "-> Will run only on data." << std::endl;
  }


  std::string outputdir = cfg.getEventYieldDir();
  system(Form("mkdir -p %s", outputdir.c_str()));


  std::string outfileName(Form("%s/trees_tmp.root", outputdir.c_str()));

  TFile* outfile = TFile::Open(outfileName.c_str(), "update");
  outfile->cd();





  if( !onlyData && !onlySignal ) { // run on MC


    std::string samplesFileName = "../samples/samples_" + cfg.mcSamples() + ".dat";
    std::cout << std::endl << std::endl;
    std::cout << "-> Loading samples from file: " << samplesFileName << std::endl;


    std::vector<ZGSample> fSamples = ZGSample::loadSamples(samplesFileName);
    //std::vector<ZGSample> fSamples = ZGSample::loadSamples(samplesFileName, 100, 999);
    if( fSamples.size()==0 ) {
      std::cout << "There must be an error: samples is empty!" << std::endl;
      exit(120);
    }


    //addTreeToFile( outfile, "zg", fSamples, cfg);
    addTreeToFile( outfile, "zg", fSamples, cfg, 851, 852 );
    addTreeToFile( outfile, "dy", fSamples, cfg, 700, 710 );
    //addTreeToFile( outfile, "top", fSamples, cfg, 300, 499 ); // irrelevant
    

    std::cout << "-> Done looping on MC samples." << std::endl;

    
  } // if MC samples



  // load signal samples, if any
  if( cfg.mcSamples()!="" && cfg.additionalStuff()!="noSignals" && !noSignals && !onlyData ) {

    std::string samplesFileName = "../samples/samples_" + cfg.mcSamples() + ".dat";
    std::cout << std::endl << std::endl;
    std::cout << "-> Loading signal samples from file: " << samplesFileName << std::endl;

    std::vector<ZGSample> fSamples = ZGSample::loadSamples(samplesFileName, 1000); // only signal (id>=1000)


    if( fSamples.size()==0 ) {

      std::cout << "No signal samples found, skipping." << std::endl;

    } else {

      for( unsigned i=0; i<fSamples.size(); ++i ) 
        addTreeToFile( outfile, fSamples[i], cfg );
    
    } // if samples != 0

  } // if mc samples

  

  if( !(cfg.dummyAnalysis()) && cfg.dataSamples()!="" && !onlyMC  && !onlySignal ) {


    std::string samplesFile_data = "../samples/samples_" + cfg.dataSamples() + ".dat";

    std::cout << std::endl << std::endl;
    std::cout << "-> Loading data from file: " << samplesFile_data << std::endl;

    std::vector<ZGSample> samples_data = ZGSample::loadSamples(samplesFile_data);
    if( samples_data.size()==0 ) {
      std::cout << "There must be an error: samples_data is empty!" << std::endl;
      exit(1209);
    }


    addTreeToFile( outfile, "data" , samples_data, cfg );
    
    std::cout << "-> Done looping on data." << std::endl;

  }


  outfile->Close();

  std::string finalFileName = outputdir + "/trees.root";
  system( Form("cp %s %s", outfileName.c_str(), finalFileName.c_str()) );

  std::cout << "-> Wrote trees to file: " << finalFileName << std::endl;



  return 0;

}




void addTreeToFile( TFile* file, ZGSample sample, const ZGConfig& cfg ) {

  std::vector<ZGSample> vec;
  vec.push_back( sample );
  return addTreeToFile( file, sample.name, vec, cfg, sample.id, sample.id );

}


void addTreeToFile( TFile* file, const std::string& treeName, std::vector<ZGSample> samples, const ZGConfig& cfg, int idMin, int idMax ) {


  if( idMin>=0 && idMax<0 ) idMax=idMin;

  TChain* tree = new TChain("mt2");


  for( unsigned i=0; i<samples.size(); ++i ) {
    if( idMin>0 ) {
      int thisId = samples[i].id;
      if( thisId<idMin ) continue;
      if( thisId>idMax ) continue;
    }
    std::string fileName(Form("%s/mt2", samples[i].file.c_str()) );
    tree->Add( fileName.c_str() );
    std::cout << "-> Added: " << fileName << std::endl;
  }
    


  //TFile* puFile_data = TFile::Open("puData.root");
  //TH1D* h1_nVert_data = (TH1D*)puFile_data->Get("nVert");

  //TFile* puFile_mc = TFile::Open("puMC.root");
  //TH1D* h1_nVert_mc = (TH1D*)puFile_mc->Get("nVert");

  // for pu reweighting:
  TFile* puFile = TFile::Open("DoubleEG_cert_210.root");
  TH1D* h1_pu = (TH1D*)puFile->Get("pileup");


  ZGTree myTree;
  myTree.loadGenStuff = false;
  myTree.Init(tree);


  // compute denominators of efficiency for PDF syst
  std::vector<float> pdf_denom;
  std::vector<float> pdf_num;

  file->cd();
  gDirectory->Delete(Form("%s;*", treeName.c_str()));
  gDirectory->Delete(Form("pdf_%s;*", treeName.c_str()));
  TTree* outTree = new TTree( treeName.c_str(), "" );

  int run;
  outTree->Branch( "run", &run, "run/I");
  int lumi;
  outTree->Branch( "lumi", &lumi, "lumi/I");
  UInt_t event;
  outTree->Branch( "event", &event, "event/i");
  float weight;
  outTree->Branch( "weight", &weight, "weight/F");
  float puWeight;
  outTree->Branch( "puWeight", &puWeight, "puWeight/F");
  int id;
  outTree->Branch( "id", &id, "id/I");
  int leptType;
  outTree->Branch( "leptType", &leptType, "leptType/I");
  float met;
  outTree->Branch( "met", &met, "met/F" );
  int nVert;
  outTree->Branch( "nVert", &nVert, "nVert/I" );
  int nGamma;
  outTree->Branch( "nGamma", &nGamma, "nGamma/I" );
  bool isGolden;
  outTree->Branch( "isGolden", &isGolden, "isGolden/O");
  bool isSilver;
  outTree->Branch( "isSilver", &isSilver, "isSilver/O");

  bool passHLT_old;
  outTree->Branch( "passHLT_old", &passHLT_old, "passHLT_old/O" );
  bool passHLT;
  outTree->Branch( "passHLT", &passHLT, "passHLT/O" );
  bool HLT_Photon165;
  outTree->Branch( "HLT_Photon165", &HLT_Photon165, "HLT_Photon165/O" );
  bool HLT_DoubleEle;
  outTree->Branch( "HLT_DoubleEle", &HLT_DoubleEle, "HLT_DoubleEle/O" );
  bool HLT_DoubleEle33;
  outTree->Branch( "HLT_DoubleEle33", &HLT_DoubleEle33, "HLT_DoubleEle33/O" );
  bool HLT_DoubleMu;
  outTree->Branch( "HLT_DoubleMu", &HLT_DoubleMu, "HLT_DoubleMu/O" );
  bool HLT_Mu30_TkMu11;
  outTree->Branch( "HLT_Mu30_TkMu11", &HLT_Mu30_TkMu11, "HLT_Mu30_TkMu11/O" );
  bool HLT_SingleMu;
  outTree->Branch( "HLT_SingleMu", &HLT_SingleMu, "HLT_SingleMu/O" );
  bool passStandardIso;
  outTree->Branch( "passStandardIso", &passStandardIso, "passStandardIso/O" );

  float weight_scale;
  outTree->Branch( "weight_scale", &weight_scale, "weight_scale/F");
  float weight_pdf;
  outTree->Branch( "weight_pdf", &weight_pdf, "weight_pdf/F");
  float weight_lep;
  outTree->Branch( "weight_lep", &weight_lep, "weight_lep/F");

  float lept0_pt;
  outTree->Branch( "lept0_pt", &lept0_pt, "lept0_pt/F" );
  float lept0_eta;
  outTree->Branch( "lept0_eta", &lept0_eta, "lept0_eta/F" );
  float lept0_phi;
  outTree->Branch( "lept0_phi", &lept0_phi, "lept0_phi/F" );
  float lept0_mass;
  outTree->Branch( "lept0_mass", &lept0_mass, "lept0_mass/F" );
  float lept0_miniRelIso;
  outTree->Branch( "lept0_miniRelIso", &lept0_miniRelIso, "lept0_miniRelIso/F" );
  int lept0_pdgId;
  outTree->Branch( "lept0_pdgId", &lept0_pdgId, "lept0_pdgId/I" );

  float lept1_pt;
  outTree->Branch( "lept1_pt", &lept1_pt, "lept1_pt/F" );
  float lept1_eta;
  outTree->Branch( "lept1_eta", &lept1_eta, "lept1_eta/F" );
  float lept1_phi;
  outTree->Branch( "lept1_phi", &lept1_phi, "lept1_phi/F" );
  float lept1_mass;
  outTree->Branch( "lept1_mass", &lept1_mass, "lept1_mass/F" );
  float lept1_miniRelIso;
  outTree->Branch( "lept1_miniRelIso", &lept1_miniRelIso, "lept1_miniRelIso/F" );
  int lept1_pdgId;
  outTree->Branch( "lept1_pdgId", &lept1_pdgId, "lept1_pdgId/I" );

  float deltaR_lept;
  outTree->Branch( "deltaR_lept", &deltaR_lept, "deltaR_lept/F" );

  float gamma_pt;
  outTree->Branch( "gamma_pt", &gamma_pt, "gamma_pt/F" );
  float gamma_eta;
  outTree->Branch( "gamma_eta", &gamma_eta, "gamma_eta/F" );
  float gamma_phi;
  outTree->Branch( "gamma_phi", &gamma_phi, "gamma_phi/F" );
  float gamma_mass;
  outTree->Branch( "gamma_mass", &gamma_mass, "gamma_mass/F" );
  float gamma_iso;
  outTree->Branch( "gamma_iso", &gamma_iso, "gamma_iso/F" );

  float z_pt;
  outTree->Branch( "z_pt", &z_pt, "z_pt/F" );
  float z_eta;
  outTree->Branch( "z_eta", &z_eta, "z_eta/F" );
  float z_phi;
  outTree->Branch( "z_phi", &z_phi, "z_phi/F" );
  float z_mass;
  outTree->Branch( "z_mass", &z_mass, "z_mass/F" );

  float boss_pt;
  outTree->Branch( "boss_pt", &boss_pt, "boss_pt/F" );
  float boss_eta;
  outTree->Branch( "boss_eta", &boss_eta, "boss_eta/F" );
  float boss_phi;
  outTree->Branch( "boss_phi", &boss_phi, "boss_phi/F" );
  float boss_mass;
  outTree->Branch( "boss_mass", &boss_mass, "boss_mass/F" );

  float boss2_pt;
  outTree->Branch( "boss2_pt", &boss2_pt, "boss2_pt/F" );
  float boss2_eta;
  outTree->Branch( "boss2_eta", &boss2_eta, "boss2_eta/F" );
  float boss2_phi;
  outTree->Branch( "boss2_phi", &boss2_phi, "boss2_phi/F" );
  float boss2_mass;
  outTree->Branch( "boss2_mass", &boss2_mass, "boss2_mass/F" );

  float qgamma_mass;
  outTree->Branch( "qgamma_mass", &qgamma_mass, "qgamma_mass/F" );
  float qZ_mass;
  outTree->Branch( "qZ_mass", &qZ_mass, "qZ_mass/F" );





  //// for muon rochester corrections
  //rochcor2015 *rmcor = new rochcor2015();

  // for photon energy scale/smearings (76X only)
  EnergyScaleCorrection_class egcor("Golden10June_plus_DCS");
  egcor.doSmearings=true;
  egcor.doScale=false;

 
  int nentries = tree->GetEntries();
  
  for( int iEntry=0; iEntry<nentries; ++iEntry ) {
    
    if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;
    
    myTree.GetEntry(iEntry);

    run   = myTree.run;
    lumi  = myTree.lumi;
    event = myTree.evt;
    id    = myTree.evt_id;

    if( event == DEBUG_EVENT ) {
      std::cout << "++++ STARTING DEBUG COUT FOR: " << std::endl;
      std::cout << "    Run: " << run << std::endl;
      std::cout << "    LS : " << lumi << std::endl;
      std::cout << "    Event : " << event << std::endl;
    }

    bool doSystForThisSample =  doSyst && (id==851 || id==4301 || id==4302);

    if( iEntry==0 && doSystForThisSample ) { // allocate all of the PDF stuff

      for( int i=0; i<myTree.nLHEweight; ++i ) {
        pdf_num.push_back(0.);
        std::cout << "[PDF Systematics] Allocating stuff for LHE weight: " << myTree.LHEweight_id[i] << " (" << i << "/" << myTree.nLHEweight << ")" << std::endl;
        bool goodIndex = (myTree.LHEweight_id[i]>=2000 && myTree.LHEweight_id[i]<=3000);
        if( goodIndex ) {
          TH1D* h1_tmp = new TH1D("pdf_tmp", "", 100, 0., 10000. );
          h1_tmp->Sumw2();
          tree->Project( "pdf_tmp", "met_pt", Form("LHEweight_wgt[%d]", i) );
          pdf_denom.push_back( h1_tmp->Integral() );
          delete h1_tmp;
        } else {
          pdf_denom.push_back( 1. );
        }
      }

    } // if first entry


    // remove overlap from DY:
    if( id>=700 && id<710 ) {
      if( myTree.ngamma>0 && myTree.gamma_mcMatchId[0]==22 ) continue;
      //if( myTree.ngamma>0 ) {
      //  bool isFake = myTree.gamma_mcMatchId[0]!=22;
      //  bool okFromDY = isFake || (!isFake && myTree.gamma_drMinParton[0]<0.05);
      //  if( !okFromDY ) continue;
      //}
    }

    if( myTree.nVert==0 ) continue;
    nVert = myTree.nVert;

    
    // filters 
    if( myTree.isData ) {
      if( !myTree.passFilters() ) continue;
    }

    isGolden = myTree.isGolden;
    isSilver = myTree.isSilver;

    HLT_Photon165   = myTree.HLT_Photon165_HE10;
    HLT_DoubleEle   = myTree.HLT_DoubleEl;
    HLT_DoubleMu    = myTree.HLT_DoubleMu;
    HLT_Mu30_TkMu11 = myTree.HLT_Mu30_TkMu11;
    HLT_DoubleEle33 = myTree.HLT_DoubleEle33;
    HLT_SingleMu    = myTree.HLT_SingleMu;

    // hlt on data:
    if( myTree.isData ) {
      if( !myTree.isGolden ) continue;
      //if( !myTree.isSilver ) continue;
      passHLT = false;
      if( id==5 ) passHLT = myTree.HLT_DoubleMu || myTree.HLT_Mu30_TkMu11; //DoubleMu PD
      if( id==8 ) passHLT = myTree.HLT_SingleMu && !myTree.HLT_DoubleMu && !myTree.HLT_Mu30_TkMu11; //SingleMu PD
      if( id==4 ) passHLT = (myTree.HLT_DoubleEl || myTree.HLT_DoubleEle33) && !myTree.HLT_DoubleMu && !myTree.HLT_Mu30_TkMu11 && !myTree.HLT_SingleMu; //DoubleEG PD
      //if( id==7 ) passHLT = myTree.HLT_Photon165_HE10 && !myTree.HLT_DoubleEl && !myTree.HLT_DoubleMu; //SinglePhoton PD
      //if( id==5 ) passHLT = myTree.HLT_DoubleMu; //DoubleMu PD
      //if( id==4 ) passHLT = myTree.HLT_DoubleEl && !myTree.HLT_DoubleMu && !myTree.HLT_SingleMu; //DoubleEG PD
      //if( id==9 ) passHLT = myTree.HLT_SingleEl && !myTree.HLT_DoubleEl && !myTree.HLT_DoubleMu && !myTree.HLT_SingleMu; //SingleElectron PD
    } else {
      // hlt on mc:
      passHLT = true; // no HLT in 80X MC
      //passHLT = ( myTree.HLT_DoubleEl || myTree.HLT_DoubleMu || myTree.HLT_Mu30_TkMu11 || myTree.HLT_DoubleEle33 );
      //passHLT = ( myTree.HLT_DoubleEl || myTree.HLT_DoubleMu || myTree.HLT_Mu30_TkMu11 || myTree.HLT_DoubleEle33 || myTree.HLT_SingleMu );
    }

    if( !passHLT && cfg.additionalStuff()!="noHLT" ) continue;

    if( event == DEBUG_EVENT )  std::cout << " -> passed HLT " << std::endl;

    if( myTree.nlep!=2 ) continue; // two leptons
    if( event == DEBUG_EVENT )  std::cout << " -> passed 2 lepton requirement " << std::endl;
    if( myTree.lep_pdgId[0] != -myTree.lep_pdgId[1] ) continue; // same flavour, opposite sign
    if( event == DEBUG_EVENT )  std::cout << " -> passed same flavor opposite charge requirement " << std::endl;

    leptType = abs(myTree.lep_pdgId[0]);
    if( leptType!=11 && leptType!=13 ) continue; // just in case
    if( event == DEBUG_EVENT )  std::cout << " -> passed electron/muon requirement " << std::endl;


      
    weight = 1.;
    puWeight = 1.;
    weight_scale = 1.;
    weight_pdf = 1.;
    weight_lep = 1.;
    
    if( !myTree.isData ) {

      weight = myTree.evt_scale1fb*cfg.lumi();
      // pu reweighting:
      //puWeight = getPUweight( h1_pu, myTree.nTrueInt );
      puWeight = myTree.puWeight;
      weight *= puWeight;

      // lepton SF:
      float leptonSF = 1.;
      if( leptType==13 ) {
        float isoSF = 0.9997;
        float idSF_loose = 0.997;
        float idSF_tight = 0.983;
        float idSF = idSF_loose*idSF_tight;
        if( myTree.lep_tightId[0]==1 && myTree.lep_tightId[1]==1 )
          idSF = 2.*idSF_loose*idSF_tight - idSF_tight*idSF_tight;

        leptonSF = isoSF*idSF;
      }
      weight *= leptonSF; 
          
      //weight *= myTree.weight_lepsf;
      //// hlt SF:
      float hltSF = (leptType==11) ? 1.0 : 0.989;
      weight *= hltSF;

      weight_lep = myTree.weight_lepsf_UP;
      weight_scale = weight; // will compute later
      weight_pdf = weight; // will compute later

    }
    


    
    TLorentzVector lept0;
    lept0.SetPtEtaPhiM( myTree.lep_pt[0], myTree.lep_eta[0], myTree.lep_phi[0], myTree.lep_mass[0] );
    TLorentzVector lept1;
    lept1.SetPtEtaPhiM( myTree.lep_pt[1], myTree.lep_eta[1], myTree.lep_phi[1], myTree.lep_mass[1] );

    if( lept0.Pt()<25. ) continue;
    if( lept1.Pt()<20. ) continue; 

    if( leptType==11 ) {
      if( myTree.lep_tightId[0]==0 || myTree.lep_tightId[1]==0 ) continue; // loose electron ID
    } else {
      //if( myTree.lep_tightId[0]==0 || myTree.lep_tightId[1]==0 ) continue; // tight muon ID on both
      if( !( (myTree.lep_tightId[0]==1 && myTree.lep_tightId[0]>=0) 
          || (myTree.lep_tightId[1]==1 && myTree.lep_tightId[1]>=0)) ) continue; // tight muon ID on at least one of the muons
    }


    // apply rochester corrections for muons:
    if( leptType==13 ) {

      // instructions taken from https://twiki.cern.ch/twiki/pub/CMS/RochcorMuon/manual_rochcor_run2.pdf
      // linked from twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/RochcorMuon
      //float qter = 1.0; 

      //if( !myTree.isData ) {
      //  rmcor->momcor_mc(lept0, myTree.lep_pdgId[0]/(abs(myTree.lep_pdgId[0])), 0, qter);
      //  rmcor->momcor_mc(lept1, myTree.lep_pdgId[1]/(abs(myTree.lep_pdgId[1])), 0, qter);
      //} else {
      //  rmcor->momcor_data(lept0, myTree.lep_pdgId[0]/(abs(myTree.lep_pdgId[0])), 0, qter);
      //  rmcor->momcor_data(lept1, myTree.lep_pdgId[1]/(abs(myTree.lep_pdgId[1])), 0, qter);
      //}

      passStandardIso = (myTree.lep_relIso04[0]<0.25 && myTree.lep_relIso04[1]<0.25);

    } else if( leptType==11 ) {

      // already applied at heppy level
      //if( !myTree.isData ) {
      //  smearEmEnergy( lept0 );
      //  smearEmEnergy( lept1 );
      //} else {
      //  applyEmEnergyScale( lept0 );
      //  applyEmEnergyScale( lept1 );
      //}

      bool passStandardIso0 = (fabs(myTree.lep_eta[0])<1.479) ? myTree.lep_relIso03[0]<0.0893 : myTree.lep_relIso03[0]<0.121;
      bool passStandardIso1 = (fabs(myTree.lep_eta[1])<1.479) ? myTree.lep_relIso03[1]<0.0893 : myTree.lep_relIso03[1]<0.121;
      
      passStandardIso = passStandardIso0 && passStandardIso1;

    }

    lept0_miniRelIso = myTree.lep_miniRelIso[0];
    lept1_miniRelIso = myTree.lep_miniRelIso[1];

    lept0_pdgId = myTree.lep_pdgId[0];
    lept1_pdgId = myTree.lep_pdgId[1];

    TLorentzVector photon;
    TLorentzVector photon2;
    bool foundSecondPhoton = false;


    if( cfg.selection()!="veryloose" ) {

      if( myTree.ngamma==0 ) continue; // photon
      if( event == DEBUG_EVENT )  std::cout << " -> passed ngamma requirement " << std::endl;

      photon  = selectPhoton(cfg, myTree, 0, lept0, lept1, egcor);

      if( photon.Pt()<1. ) continue;
      if( event == DEBUG_EVENT )  std::cout << " -> passed photon requirement " << std::endl;

      photon2 = selectPhoton(cfg, myTree, 1, lept0, lept1, egcor);

      foundSecondPhoton = photon2.Pt()>10.;

    }


    TLorentzVector zBoson = lept0+lept1;
    if( zBoson.M()<50. ) continue;
    if( cfg.selection()!="presel" && zBoson.M()>130. ) continue;
    if( event == DEBUG_EVENT )  std::cout << " -> passed M(ll) cut " << std::endl;

    TLorentzVector boss = zBoson + photon;
    TLorentzVector boss2;
    if( foundSecondPhoton ) 
      boss2 = zBoson + photon + photon2;
    else
      boss2.SetPtEtaPhiM( 0.1, 0., 0., 0.);


    lept0_pt   = lept0.Pt();
    lept0_eta  = lept0.Eta();
    lept0_phi  = lept0.Phi();
    lept0_mass = lept0.M();

    lept1_pt   = lept1.Pt();
    lept1_eta  = lept1.Eta();
    lept1_phi  = lept1.Phi();
    lept1_mass = lept1.M();

    deltaR_lept = lept0.DeltaR(lept1);

    nGamma = myTree.ngamma;

    if( photon.Pt()>0. ) {

      gamma_pt   = photon.Pt();
      gamma_eta  = photon.Eta();
      gamma_phi  = photon.Phi();
      gamma_mass = photon.M();
      gamma_iso = myTree.gamma_chHadIso[0];

    } else {

      gamma_pt   = 0.;
      gamma_eta  = 0.;
      gamma_phi  = 0.;
      gamma_mass = 0.;
      gamma_iso = -1.;

    }

    z_pt   = zBoson.Pt();
    z_eta  = zBoson.Eta();
    z_phi  = zBoson.Phi();
    z_mass = zBoson.M();

    boss_pt   = boss.Pt();
    boss_eta  = boss.Eta();
    boss_phi  = boss.Phi();
    boss_mass = boss.M();

    boss2_pt   = boss2.Pt();
    boss2_eta  = boss2.Eta();
    boss2_phi  = boss2.Phi();
    boss2_mass = boss2.M();

    met = myTree.met_pt;

    if( cfg.selection()=="v0" )
      if( gamma_pt/boss_mass < 40./150. ) continue;

    if( cfg.selection()=="Qv0" ) {
      if( myTree.njet==0 ) continue;
      if( fabs(myTree.jet_eta[0])>2.5 ) continue;
      if( myTree.jet_pt[0]<30. ) continue;
    }

    if( event == DEBUG_EVENT )  std::cout << " -> passed additional cuts" << std::endl;


    qgamma_mass = -1.;
    qZ_mass = -1.;

    if( myTree.njet>0 && myTree.jet_pt[0]>30. ) {

      TLorentzVector jet;
      jet.SetPtEtaPhiM( myTree.jet_pt[0], myTree.jet_eta[0], myTree.jet_phi[0], myTree.jet_mass[0] );

      TLorentzVector qgamma = jet+photon;
      TLorentzVector qZ     = jet+zBoson;

      qgamma_mass = qgamma.M();
      qZ_mass     = qZ.M();

    }


    if( DATABLINDING && myTree.isData && boss_mass>500. ) continue;


    if( doSystForThisSample ) { // systematic uncertainties

      for( int i=0; i<myTree.nLHEweight; ++i ) {
        bool goodIndex = (myTree.LHEweight_id[i]>=2000 && myTree.LHEweight_id[i]<=3000);
        if( goodIndex )
          pdf_num[i] += myTree.LHEweight_wgt[i];
      } //for pdf sets 

    }


    //if( id==851 && doSyst ) { // systematic uncertainties

    //  // first scale
    //  float ref = myTree.LHEweight_original;

    //  float maxScaleDiff = 0.;

    //  TH1D* h1_pdf = new TH1D("pdf", "", 1000, -3000., 3000.);

    //  for( int i=0; i<myTree.nLHEweight; ++i ) {

    //    if( myTree.LHEweight_id[i]>=1000 && myTree.LHEweight_id[i]<1010 ) {
    //      float thisScaleDiff = fabs( (myTree.LHEweight_wgt[i]-ref)/ref );
    //      if( thisScaleDiff>maxScaleDiff) 
    //        maxScaleDiff = thisScaleDiff+1.;
    //    }

    //    if( myTree.LHEweight_id[i]>=2000 ) {
    //      if( myTree.LHEweight_wgt[i] > h1_pdf->GetXaxis()->GetXmax() || myTree.LHEweight_wgt[i] < h1_pdf->GetXaxis()->GetXmin() )
    //        std::cout << "WARNING!! PDF weight out of bounds: " << myTree.LHEweight_wgt[i] << std::endl;
    //      h1_pdf->Fill( myTree.LHEweight_wgt[i]/ref );
    //    }

    //  } // for lhe weights

    //  float meanPdfWgt = h1_pdf->GetMean();
    //  float rmsPdfWgt = h1_pdf->GetRMS();
    //  weight_pdf *= (1.+rmsPdfWgt/meanPdfWgt);

    //  weight_scale *= maxScaleDiff;

    //  delete h1_pdf;

    //} // if correct sample

    if( event == DEBUG_EVENT )  std::cout << " -> filling tree" << std::endl;
    outTree->Fill();

    
  } // for entries
    

  TH1D* h1_pdf = new TH1D( Form("pdf_%s", treeName.c_str()), "", 2000, 0., 1. );
  for( unsigned i=0; i<pdf_num.size(); ++i ) {
    if( pdf_num[i]!=0. )
      h1_pdf->Fill( pdf_num[i]/pdf_denom[i] );
  }

  if( h1_pdf->GetMean()>0. )
    h1_pdf->Write();
 
  outTree->Write();

  puFile->Close();

}



void smearEmEnergy( TLorentzVector& p ) {

  float smearEBlowEta     = 0.013898;
  float smearEBhighEta    = 0.0189895;
  float smearEElowEta     = 0.027686;
  float smearEEhighEta    = 0.031312;
  float theGaussMean      = 1.;
  float pt = p.Pt();
  float eta = p.Eta();
  float phi = p.Phi();
  float mass = p.M();
  float theSmear = 0.;
  if      (fabs(eta)<1.                   )  theSmear = smearEBlowEta;
  else if (fabs(eta)>=1.  && fabs(eta)<1.5)  theSmear = smearEBhighEta;
  else if (fabs(eta)>=1.5 && fabs(eta)<2. )  theSmear = smearEElowEta;
  else if (fabs(eta)>=2.  && fabs(eta)<2.5)  theSmear = smearEEhighEta;
  float fromGauss = myRandom_.Gaus(theGaussMean,theSmear);
  p.SetPtEtaPhiM( fromGauss*pt, eta, phi, mass ); // keep mass and direction same

}


void applyEmEnergyScale( TLorentzVector& p ) {

  float scaleEBlowEta  = 2./(0.99544 + 0.99882);
  float scaleEBhighEta = 2./(0.99662 + 1.0065);
  float scaleEElowEta  = 2./(0.98633 + 0.99536);
  float scaleEEhighEta = 2./(0.97859 + 0.98567);
  float pt = p.Pt();
  float eta = p.Eta();
  float phi = p.Phi();
  float mass = p.M();
  float theScale = 1.;
  if      (fabs(eta)<1.                     ) theScale = scaleEBlowEta;
  else if (fabs(eta)>=1.  && fabs(eta)<1.5  ) theScale = scaleEBhighEta;
  else if (fabs(eta)>=1.5 && fabs(eta)<2.   ) theScale = scaleEElowEta;
  else if (fabs(eta)>=2.  && fabs(eta)<2.5  ) theScale = scaleEEhighEta;
  p.SetPtEtaPhiM( theScale*pt, eta, phi, mass ); // keep mass and direction same

}


TLorentzVector selectPhoton( const ZGConfig& cfg, const ZGTree& myTree, int index, const TLorentzVector& lept0, const TLorentzVector& lept1, EnergyScaleCorrection_class egcor ) {

  TLorentzVector photon;
  if( myTree.ngamma<=index ) {
    photon.SetPtEtaPhiM( 0.01, 0., 0., 0. );
    return photon;
  }
  photon.SetPtEtaPhiM( myTree.gamma_pt[index], myTree.gamma_eta[index], myTree.gamma_phi[index], myTree.gamma_mass[index] );

  //// photon energy corrections/smearing 
  //if( !myTree.isData ) {
  //  float smearSigma =  egcor.getSmearingSigma(myTree.run, fabs(photon.Eta())<1.479, myTree.gamma_r9[0], photon.Eta(), photon.Pt(), 0., 0.);
  //  std::cout << "eta: " << photon.Eta() << " pt: " << photon.Pt() << " r9: " << myTree.gamma_r9[0] << " smear: " << smearSigma << std::endl;
  //  exit(11);
  //  //smearEmEnergy( photon ); // 74X smearings
  //} else {
  //  //applyEmEnergyScale( photon );
  //  float cor = egcor.ScaleCorrection(myTree.run, fabs(photon.Eta())<1.479, myTree.gamma_r9[0], photon.Eta(), photon.Pt());
  //  photon.SetPtEtaPhiM( photon.Pt()*cor, photon.Eta(), photon.Phi(), photon.M() );
  //}

  bool goodPhoton = true;
  if( photon.Pt()<40. ) goodPhoton=false;
  if( fabs(photon.Eta())>1.4442 && fabs(photon.Eta())<1.566 ) goodPhoton=false;
  if( fabs(photon.Eta())>2.5 ) goodPhoton=false;
  if( myTree.gamma_idCutBased[index]==0 ) goodPhoton=false;
  if( myTree.gamma_chHadIso[index]>2.5 ) goodPhoton=false;
  if( fabs(myTree.gamma_eta[index])<1.44 ) {
    if( myTree.gamma_sigmaIetaIeta[index]>0.0102 ) goodPhoton=false;
  } else {
    if( myTree.gamma_sigmaIetaIeta[index]>0.0274 ) goodPhoton=false;
  }
  float deltaR_thresh = 0.4;
  if( photon.DeltaR(lept0)<deltaR_thresh || photon.DeltaR(lept1)<deltaR_thresh ) goodPhoton=false;


  if( !goodPhoton )
    photon.SetPtEtaPhiM( 0.01, 0., 0., 0. );

  return photon;

}


float getPUweight( TH1D* h1_pu, float nTrueInt ) {

  int bin = h1_pu->FindBin(nTrueInt);
  float w = h1_pu->GetBinContent(bin)/h1_pu->Integral();

  return w;

}
