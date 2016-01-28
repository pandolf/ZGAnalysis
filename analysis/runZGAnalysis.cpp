#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>


#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"


#include "interface/ZGSample.h"
#include "interface/ZGConfig.h"
#include "interface/ZGCommonTools.h"



#define zg_cxx
#include "interface/zg.h"


#define DATABLINDING true




void addTreeToFile( TFile* file, const std::string& treeName, std::vector<ZGSample> samples, const ZGConfig& cfg, int idMin=-1, int idMax=-1 );




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
  if( argc > 2 ) {

    std::string dataMC(argv[2]);
    if( dataMC=="data" ) onlyData = true;
    else if( dataMC=="MC" || dataMC=="mc" ) onlyMC = true;
    else if( dataMC=="signal" ) onlySignal = true;
    else {
      std::cout << "-> You passed a second argument that isn't 'data', nor 'MC', nor 'signal', so I don't know what to do about it." << std::endl;
    }

  } else {

    std::cout << "-> Will run on both data and MC." << std::endl;

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

    std::vector<ZGSample> fSamples = ZGSample::loadSamples(samplesFileName, 101, 999); // not interested in signal here (see later)
    if( fSamples.size()==0 ) {
      std::cout << "There must be an error: samples is empty!" << std::endl;
      exit(120);
    }


    addTreeToFile( outfile, "zg", fSamples, cfg, 851 );
    addTreeToFile( outfile, "dy", fSamples, cfg, 701 );
    //addTreeToFile( outfile, "top", fSamples, cfg, 300, 499 ); // irrelevant
    

    std::cout << "-> Done looping on MC samples." << std::endl;

    
  } // if MC samples



//// load signal samples, if any
//std::vector< MT2Analysis< MT2EstimateSigSyst>* > signals;
//if( cfg.mcSamples()!="" && cfg.additionalStuff()!="noSignals" && !onlyData ) {

//  std::string samplesFileName = "../samples/samples_" + cfg.mcSamples() + ".dat";
//  std::cout << std::endl << std::endl;
//  std::cout << "-> Loading signal samples from file: " << samplesFileName << std::endl;

//  std::vector<ZGSample> fSamples = ZGSample::loadSamples(samplesFileName, 1000); // only signal (id>=1000)


//  if( fSamples.size()==0 ) {

//    std::cout << "No signal samples found, skipping." << std::endl;

//  } else {
//  
//    for( unsigned i=0; i<fSamples.size(); ++i ) 
//      signals.push_back( computeSigYield<MT2EstimateSigSyst>( fSamples[i], cfg ) );
//  
//    std::cout << "     merging T1bbbb full scan..." << std::endl;
//    MT2Analysis<MT2EstimateSigSyst>* EventYield_T1bbbb   = mergeYields<MT2EstimateSigSyst>( signals, cfg.regionsSet(), "SMS_T1bbbb_fullScan", 1020, 1020 );
//    std::cout << "-> Done merging." << std::endl;

//    signals.push_back( EventYield_T1bbbb );

//  } // if samples != 0

//} // if mc samples

//else if ( cfg.sigSamples()!="" && !onlyData ) {

//  std::string samplesFileName = "../samples/samples_" + cfg.sigSamples() + ".dat";
//  std::cout << std::endl << std::endl;
//  std::cout << "-> Loading signal samples from file: " << samplesFileName << std::endl;

//  std::vector<ZGSample> fSamples = ZGSample::loadSamples(samplesFileName, 1000); // only signal (id>=1000)

//  if( fSamples.size()==0 ) {

//    std::cout << "No signal samples found, skipping." << std::endl;

//  } else {

//    for( unsigned i=0; i<fSamples.size(); ++i )
//      signals.push_back( computeSigYield<MT2EstimateSigSyst>( fSamples[i], cfg ) );

//    std::cout << "     merging T1bbbb full scan..." << std::endl;
//    MT2Analysis<MT2EstimateSigSyst>* EventYield_T1bbbb   = mergeYields<MT2EstimateSigSyst>( signals, cfg.regionsSet(), "SMS_T1bbbb_fullScan", 1020, 1020 );
//    std::cout << "-> Done merging." << std::endl;

//    signals.push_back( EventYield_T1bbbb );

//  } // if samples != 0
//  
//} // if sig samples
//

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


//if( yields.size()==0 && signals.size()==0 ) {
//  std::cout << "-> Didn't end up with a single yield... something's wrong." << std::endl;
//  exit(87);
//}


//// save MT2Analyses:
//if( yields.size()>0 ){
//  yields[0]->writeToFile(outputdir + "/analyses.root");
//for( unsigned i=1; i<yields.size(); ++i )
//  yields[i]->writeToFile(outputdir + "/analyses.root");
//for( unsigned i=0; i<signals.size(); ++i )
//  signals[i]->writeToFile(outputdir + "/analyses.root");
//}
//else if( signals.size()>0 ){
//  signals[0]->writeToFile(outputdir + "/analyses.root");
//  for( unsigned i=1; i<signals.size(); ++i )
//    signals[i]->writeToFile(outputdir + "/analyses.root");
//}
//cfg.saveAs(outputdir + "/config.txt");


  outfile->Close();

  // move in position
  std::string finalFileName = outputdir + "/trees.root";
  system( Form("cp %s %s", outfileName.c_str(), finalFileName.c_str()) );

  std::cout << "-> Wrote trees to file: " << finalFileName << std::endl;



  return 0;

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
    tree->Add( Form("%s/mt2", samples[i].file.c_str()) );
  }
    


  TFile* puFile_data = TFile::Open("puData.root");
  TH1D* h1_nVert_data = (TH1D*)puFile_data->Get("nVert");

  TFile* puFile_mc = TFile::Open("puMC.root");
  TH1D* h1_nVert_mc = (TH1D*)puFile_mc->Get("nVert");


  ZGTree myTree;
  myTree.loadGenStuff = false;
  myTree.Init(tree);



  file->cd();
  gDirectory->Delete(Form("%s;*", treeName.c_str()));
  TTree* outTree = new TTree( treeName.c_str(), "" );

  int run;
  outTree->Branch( "run", &run, "run/I");
  int lumi;
  outTree->Branch( "lumi", &lumi, "lumi/I");
  unsigned event;
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

  float lept0_pt;
  outTree->Branch( "lept0_pt", &lept0_pt, "lept0_pt/F" );
  float lept0_eta;
  outTree->Branch( "lept0_eta", &lept0_eta, "lept0_eta/F" );
  float lept0_phi;
  outTree->Branch( "lept0_phi", &lept0_phi, "lept0_phi/F" );

  float lept1_pt;
  outTree->Branch( "lept1_pt", &lept1_pt, "lept1_pt/F" );
  float lept1_eta;
  outTree->Branch( "lept1_eta", &lept1_eta, "lept1_eta/F" );
  float lept1_phi;
  outTree->Branch( "lept1_phi", &lept1_phi, "lept1_phi/F" );

  float gamma_pt;
  outTree->Branch( "gamma_pt", &gamma_pt, "gamma_pt/F" );
  float gamma_eta;
  outTree->Branch( "gamma_eta", &gamma_eta, "gamma_eta/F" );
  float gamma_phi;
  outTree->Branch( "gamma_phi", &gamma_phi, "gamma_phi/F" );
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

 
  int nentries = tree->GetEntries();
  
  for( int iEntry=0; iEntry<nentries; ++iEntry ) {
    
    if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;
    
    myTree.GetEntry(iEntry);

    run   = myTree.run;
    lumi  = myTree.lumi;
    event = myTree.evt;
    id    = myTree.evt_id;

    // remove overlap
    if( id==701 && myTree.ngamma>0 && myTree.gamma_mcMatchId[0]==22 ) continue;

    if( myTree.nVert==0 ) continue;
    nVert = myTree.nVert;

    
    // filters 
    if( myTree.isData ) {
      if( !myTree.passFilters() ) continue;
    }

    // hlt
    if(  myTree.isData && !( ( id==5  && myTree.HLT_DoubleMu ) || ( id==4 && myTree.HLT_DoubleEl ) || ( id==7 && !myTree.HLT_DoubleEl && !myTree.HLT_DoubleMu && myTree.HLT_Photon165_HE10 ) )  ) continue;
    //if(  myTree.isData && !( ( id==5  && myTree.HLT_DoubleMu ) || ( id==4 && (myTree.HLT_DoubleEl || myTree.HLT_SingleEl) ) || ( id==7 && !myTree.HLT_DoubleEl && !myTree.HLT_DoubleMu && myTree.HLT_Photon165_HE10 && !myTree.HLT_SingleEl ) )  ) continue;
    //if(  myTree.isData && !( ( id==5  && myTree.HLT_DoubleMu ) || ( id==4 ) || ( id==7 && !myTree.HLT_DoubleEl && !myTree.HLT_DoubleMu && myTree.HLT_Photon165_HE10 ) )  ) continue;

    //if( !myTree.isData && !( myTree.HLT_DoubleEl || myTree.HLT_DoubleMu )) continue;
    if( !myTree.isData && !( myTree.HLT_DoubleEl || myTree.HLT_DoubleMu || myTree.HLT_Photon165_HE10)) continue;
    //if( !myTree.isData && !( myTree.HLT_DoubleEl || myTree.HLT_DoubleMu || myTree.HLT_Photon165_HE10 || myTree.HLT_SingleEl)) continue;
      
    weight = 1.;
    puWeight = 1.;
    // pu reweighting:
    if( !myTree.isData ) {
      puWeight = ZGCommonTools::getPUweight( nVert, h1_nVert_data, h1_nVert_mc );
      weight = myTree.evt_scale1fb*puWeight;
    }
    

    if( myTree.nlep!=2 ) continue; // two leptons
    if( myTree.lep_pdgId[0] != -myTree.lep_pdgId[1] ) continue; // same flavour, opposite sign

    leptType = abs(myTree.lep_pdgId[0]);
    if( leptType!=11 && leptType!=13 ) continue; // just in case

    
    TLorentzVector lept0;
    lept0.SetPtEtaPhiM( myTree.lep_pt[0], myTree.lep_eta[0], myTree.lep_phi[0], myTree.lep_mass[0] );
    TLorentzVector lept1;
    lept1.SetPtEtaPhiM( myTree.lep_pt[1], myTree.lep_eta[1], myTree.lep_phi[1], myTree.lep_mass[1] );

    if( lept0.Pt()<25. ) continue;
    if( lept1.Pt()<20. ) continue; 

    if( leptType==11 ) {
      if( myTree.lep_tightId[0]==0 || myTree.lep_tightId[1]==0 ) continue; // loose electron ID
    }

    TLorentzVector photon;


    if( cfg.selection()!="veryloose" ) {

      if( myTree.ngamma==0 ) continue; // photon
      photon.SetPtEtaPhiM( myTree.gamma_pt[0], myTree.gamma_eta[0], myTree.gamma_phi[0], myTree.gamma_mass[0] );
   
      if( photon.Pt()<40. ) continue;
      if( fabs(photon.Eta())>1.44 && fabs(photon.Eta())<1.57 ) continue;
      if( fabs(photon.Eta())>2.5 ) continue;
      if( myTree.gamma_idCutBased[0]==0 ) continue;
      if( myTree.gamma_chHadIso[0]>2.5 ) continue;
      float deltaR_thresh = 0.4;
      if( photon.DeltaR(lept0)<deltaR_thresh || photon.DeltaR(lept1)<deltaR_thresh ) continue;

    }


    TLorentzVector zBoson = lept0+lept1;
    if( zBoson.M()<50. ) continue;

    TLorentzVector boss = zBoson + photon;


    lept0_pt  = lept0.Pt();
    lept0_eta = lept0.Eta();
    lept0_phi = lept0.Phi();

    lept1_pt  = lept1.Pt();
    lept1_eta = lept1.Eta();
    lept1_phi = lept1.Phi();

    nGamma = myTree.ngamma;

    if( photon.Pt()>0. ) {

      gamma_pt  = photon.Pt();
      gamma_eta = photon.Eta();
      gamma_phi = photon.Phi();
      gamma_iso = myTree.gamma_chHadIso[0];

    } else {

      gamma_pt  = 0.;
      gamma_eta = 0.;
      gamma_phi = 0.;
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

    met = myTree.met_pt;
    //if( met > 80. ) continue;
    if( gamma_pt/boss_mass< 40./150. ) continue;

    if( DATABLINDING && myTree.isData && boss_mass>500. ) continue;

    outTree->Fill();
    
  } // for entries
    
  
  outTree->Write();

}



