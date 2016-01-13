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



#define zg_cxx
#include "interface/zg.h"





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



  if( !onlyData && !onlySignal ) { // run on MC

    std::string outfileName(Form("%s/mc.root", outputdir.c_str()));
    TFile* outfile = TFile::Open(outfileName.c_str(), "recreate");
    outfile->cd();


    std::string samplesFileName = "../samples/samples_" + cfg.mcSamples() + ".dat";
    std::cout << std::endl << std::endl;
    std::cout << "-> Loading samples from file: " << samplesFileName << std::endl;

    std::vector<ZGSample> fSamples = ZGSample::loadSamples(samplesFileName, 101, 999); // not interested in signal here (see later)
    if( fSamples.size()==0 ) {
      std::cout << "There must be an error: samples is empty!" << std::endl;
      exit(120);
    }


    addTreeToFile( outfile, "zg" , fSamples, cfg, 851 );
    addTreeToFile( outfile, "top", fSamples, cfg, 300, 499 );
    
    outfile->Close();

    std::cout << "-> Done looping on MC samples." << std::endl;
    std::cout << "-> Saved trees in: " << outfileName << std::endl;

    
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

//if( !(cfg.dummyAnalysis()) && cfg.dataSamples()!="" && !onlyMC  && !onlySignal ) {

//  std::string samplesFile_data = "../samples/samples_" + cfg.dataSamples() + ".dat";

//  std::cout << std::endl << std::endl;
//  std::cout << "-> Loading data from file: " << samplesFile_data << std::endl;

//  //    std::vector<ZGSample> samples_data = ZGSample::loadSamples(samplesFile_data, "JetHTMHT"); //, 1, 99 );
//  //    std::vector<ZGSample> samples_data = ZGSample::loadSamples(samplesFile_data, 1, 3 );
//  std::vector<ZGSample> samples_data = ZGSample::loadSamples(samplesFile_data, -1, 0 );
//  if( samples_data.size()==0 ) {
//    std::cout << "There must be an error: samples_data is empty!" << std::endl;
//    exit(1209);
//  }

//  std::vector< MT2Analysis<MT2EstimateTree>* > EventYield_data;
//  for( unsigned i=0; i < samples_data.size(); ++i )
//    EventYield_data.push_back( computeYield<MT2EstimateTree>( samples_data[i], cfg ) );

//  MT2Analysis<MT2EstimateTree>* dataYield;
//  //dataYield = EventYield_data[0];
//  //dataYield->setName("data");
//  //dataYield   = mergeYields<MT2EstimateTree>( EventYield_data, cfg.regionsSet(), "data", 1, 3 );
//  dataYield   = mergeYields<MT2EstimateTree>( EventYield_data, cfg.regionsSet(), "data", -1, -1 );

//  yields.push_back( dataYield );

//}


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
    


  ZGTree myTree;
  myTree.loadGenStuff = false;
  myTree.Init(tree);


  file->cd();

  TTree* outTree = new TTree( treeName.c_str(), "" );

  int run;
  outTree->Branch( "run", &run, "run/I");
  int lumi;
  outTree->Branch( "lumi", &lumi, "lumi/I");
  unsigned event;
  outTree->Branch( "event", &event, "event/i");
  float weight;
  outTree->Branch( "weight", &weight, "weight/F");
  int id;
  outTree->Branch( "id", &id, "id/I");
  int leptType;
  outTree->Branch( "leptType", &leptType, "leptType/I");

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
    
    // filters
    if( myTree.isData ) {
      if( !myTree.passFilters() ) continue;
    }

    run   = myTree.run;
    lumi  = myTree.lumi;
    event = myTree.evt;
    id    = myTree.evt_id;

    // hlt
    if(  myTree.isData && !( ( id==5  && myTree.HLT_DoubleMu ) || ( id==4 && myTree.HLT_DoubleEl ) || ( id==7 && !myTree.HLT_DoubleEl && !myTree.HLT_DoubleMu && myTree.HLT_Photon165_HE10 ) )  ) continue;
    if( !myTree.isData && !( myTree.HLT_DoubleEl || myTree.HLT_DoubleMu || myTree.HLT_Photon165_HE10)) continue;
      
    weight = (myTree.isData) ? 1. : myTree.evt_scale1fb*cfg.lumi();

    if( myTree.nlep!=2 ) continue; // two leptons
    if( myTree.lep_pdgId[0] != -myTree.lep_pdgId[1] ) continue; // same flavour, opposite sign

    leptType = abs(myTree.lep_pdgId[0]);
    if( leptType!=11 && leptType!=13 ) continue; // just in case

    if( myTree.ngamma==0 ) continue; // photon
    
    TLorentzVector lept0;
    lept0.SetPtEtaPhiM( myTree.lep_pt[0], myTree.lep_eta[0], myTree.lep_phi[0], myTree.lep_mass[0] );
    TLorentzVector lept1;
    lept1.SetPtEtaPhiM( myTree.lep_pt[1], myTree.lep_eta[1], myTree.lep_phi[1], myTree.lep_mass[1] );

    if( lept0.Pt()<25. ) continue;
    if( lept1.Pt()<20. ) continue; 


    TLorentzVector photon;
    photon.SetPtEtaPhiM( myTree.gamma_pt[0], myTree.gamma_eta[0], myTree.gamma_phi[0], myTree.gamma_mass[0] );

    if( photon.Pt()<50. ) continue;
    if( fabs(photon.Eta())>1.44 && fabs(photon.Eta())<1.57 ) continue;
    if( myTree.gamma_chHadIso[0]>2.5 ) continue;
    if( photon.DeltaR(lept0)<0.5 || photon.DeltaR(lept1)<0.5 ) continue;


    TLorentzVector zBoson = lept0+lept1;
    if( zBoson.M()<50. ) continue;

    TLorentzVector boss = zBoson + photon;


    lept0_pt  = lept0.Pt();
    lept0_eta = lept0.Eta();
    lept0_phi = lept0.Phi();

    lept1_pt  = lept1.Pt();
    lept1_eta = lept1.Eta();
    lept1_phi = lept1.Phi();

    gamma_pt  = photon.Pt();
    gamma_eta = photon.Eta();
    gamma_phi = photon.Phi();

    z_pt   = zBoson.Pt();
    z_eta  = zBoson.Eta();
    z_phi  = zBoson.Phi();
    z_mass = zBoson.M();

    boss_pt   = boss.Pt();
    boss_eta  = boss.Eta();
    boss_phi  = boss.Phi();
    boss_mass = boss.M();


    outTree->Fill();
    
  } // for entries
    
  
  outTree->Write();

}

