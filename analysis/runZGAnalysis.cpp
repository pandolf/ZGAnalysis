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

bool doSyst = true;


void addTreeToFile( TFile* file, ZGSample sample, const ZGConfig& cfg );
void addTreeToFile( TFile* file, const std::string& treeName, std::vector<ZGSample> samples, const ZGConfig& cfg, int idMin=-1, int idMax=-1 );
void smearEmEnergy     ( TLorentzVector& p );
void applyEmEnergyScale( TLorentzVector& p );



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
    else if( dataMC=="mcbg" || dataMC=="mcBG" || dataMC=="MCBG" ) {
      onlyMC = true;
      noSignals = true;
    } else if( dataMC=="signal" ) onlySignal = true;
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
  bool HLT_SingleMu;
  outTree->Branch( "HLT_SingleMu", &HLT_SingleMu, "HLT_SingleMu/O" );

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

  float lept1_pt;
  outTree->Branch( "lept1_pt", &lept1_pt, "lept1_pt/F" );
  float lept1_eta;
  outTree->Branch( "lept1_eta", &lept1_eta, "lept1_eta/F" );
  float lept1_phi;
  outTree->Branch( "lept1_phi", &lept1_phi, "lept1_phi/F" );

  float deltaR_lept;
  outTree->Branch( "deltaR_lept", &deltaR_lept, "deltaR_lept/F" );

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



  // for muon rochester corrections
  rochcor2015 *rmcor = new rochcor2015();

  // for photon energy scale/smearings (76X only)
  EnergyScaleCorrection_class egcor("76X_16DecRereco_2015");
  egcor.doScale=true;
  egcor.doSmearings=false;

 
  int nentries = tree->GetEntries();
  
  for( int iEntry=0; iEntry<nentries; ++iEntry ) {
    
    if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;
    
    myTree.GetEntry(iEntry);

    run   = myTree.run;
    lumi  = myTree.lumi;
    event = myTree.evt;
    id    = myTree.evt_id;

    // remove overlap from DY:
    if( id>=700 && id<710 ) {
      if( myTree.ngamma>0 && myTree.gamma_mcMatchId[0]==22 ) continue;
      //if( myTree.ngamma==0 ) continue;
      //bool isFake = myTree.gamma_mcMatchId[0]!=22;
      //bool okFromDY = isFake || (!isFake && myTree.gamma_drMinParton[0]<0.05);
      //if( !okFromDY ) continue;
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
    HLT_DoubleEle33 = myTree.HLT_DoubleEl33;
    HLT_SingleMu    = myTree.HLT_SingleMu;

    // hlt on data:
    if( myTree.isData ) {
      //if( !myTree.isGolden ) continue;
      if( !myTree.isSilver ) continue;
      passHLT = false;
      if( id==5 ) passHLT = myTree.HLT_DoubleMu; //DoubleMu PD
      if( id==8 ) passHLT = myTree.HLT_SingleMu && !myTree.HLT_DoubleMu; //SingleMu PD
      if( id==4 ) passHLT = (myTree.HLT_DoubleEl || myTree.HLT_DoubleEl33) && !myTree.HLT_DoubleMu && !myTree.HLT_SingleMu; //DoubleEG PD
      //if( id==7 ) passHLT = myTree.HLT_Photon165_HE10 && !myTree.HLT_DoubleEl && !myTree.HLT_DoubleMu; //SinglePhoton PD
      //if( id==5 ) passHLT = myTree.HLT_DoubleMu; //DoubleMu PD
      //if( id==4 ) passHLT = myTree.HLT_DoubleEl && !myTree.HLT_DoubleMu && !myTree.HLT_SingleMu; //DoubleEG PD
      //if( id==9 ) passHLT = myTree.HLT_SingleEl && !myTree.HLT_DoubleEl && !myTree.HLT_DoubleMu && !myTree.HLT_SingleMu; //SingleElectron PD
    } else {
      // hlt on mc:
      //passHLT = ( myTree.HLT_DoubleEl || myTree.HLT_DoubleMu || myTree.HLT_DoubleEl33 );
      passHLT = ( myTree.HLT_DoubleEl || myTree.HLT_DoubleMu || myTree.HLT_DoubleEl33 || myTree.HLT_SingleMu );
    }

    if( !passHLT ) continue;

    if( myTree.nlep!=2 ) continue; // two leptons
    if( myTree.lep_pdgId[0] != -myTree.lep_pdgId[1] ) continue; // same flavour, opposite sign

    leptType = abs(myTree.lep_pdgId[0]);
    if( leptType!=11 && leptType!=13 ) continue; // just in case


      
    weight = 1.;
    puWeight = 1.;
    weight_scale = 1.;
    weight_pdf = 1.;
    weight_lep = 1.;
    
    if( !myTree.isData ) {

      weight = myTree.evt_scale1fb;
      // pu reweighting:
      puWeight = myTree.puWeight;
      //puWeight = ZGCommonTools::getPUweight( nVert, h1_nVert_data, h1_nVert_mc );
      weight *= puWeight;
      // lepton SF:
      weight *= myTree.weight_lepsf;
      //// hlt SF:
      //float hltSF = (leptType==11) ? 1.02 : 0.94;
      //weight *= hltSF;

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
      if( myTree.lep_tightId[0]==0 && myTree.lep_tightId[1]==0 ) continue; // tight muon ID on at least one of the muons
    }


    // apply rochester corrections for muons:
    if( leptType==13 ) {

      // instructions taken from https://twiki.cern.ch/twiki/pub/CMS/RochcorMuon/manual_rochcor_run2.pdf
      // linked from twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/RochcorMuon
      float qter = 1.0; 

      if( !myTree.isData ) {
        rmcor->momcor_mc(lept0, myTree.lep_pdgId[0]/(abs(myTree.lep_pdgId[0])), 0, qter);
        rmcor->momcor_mc(lept1, myTree.lep_pdgId[1]/(abs(myTree.lep_pdgId[1])), 0, qter);
      } else {
        rmcor->momcor_data(lept0, myTree.lep_pdgId[0]/(abs(myTree.lep_pdgId[0])), 0, qter);
        rmcor->momcor_data(lept1, myTree.lep_pdgId[1]/(abs(myTree.lep_pdgId[1])), 0, qter);
      }

    } else if( leptType==11 ) {

      // already applied at heppy level
      //if( !myTree.isData ) {
      //  smearEmEnergy( lept0 );
      //  smearEmEnergy( lept1 );
      //} else {
      //  applyEmEnergyScale( lept0 );
      //  applyEmEnergyScale( lept1 );
      //}
      
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
      if( fabs(myTree.gamma_eta[0])<1.44 ) {
        if( myTree.gamma_sigmaIetaIeta[0]>0.0102 ) continue;
      } else {
        if( myTree.gamma_sigmaIetaIeta[0]>0.0274 ) continue;
      }
      float deltaR_thresh = 0.4;
      if( photon.DeltaR(lept0)<deltaR_thresh || photon.DeltaR(lept1)<deltaR_thresh ) continue;

      // photon energy corrections/smearing NOT done at heppy (in 74X, will be in for 76X)
      if( !myTree.isData ) {
        if( cfg.smearing() )
          smearEmEnergy( photon );
      } else {
        applyEmEnergyScale( photon );
        //float cor = egcor.ScaleCorrection(run, fabs(photon.Eta())<1.479, myTree.gamma_r9[0], photon.Eta(), photon.Pt());
        //photon.SetPtEtaPhiM( photon.Pt()*cor, photon.Eta(), photon.Phi(), photon.M() );
      }

    }


    TLorentzVector zBoson = lept0+lept1;
    if( zBoson.M()<50. ) continue;
    if( cfg.selection()!="presel" && zBoson.M()>130. ) continue;

    TLorentzVector boss = zBoson + photon;


    lept0_pt  = lept0.Pt();
    lept0_eta = lept0.Eta();
    lept0_phi = lept0.Phi();

    lept1_pt  = lept1.Pt();
    lept1_eta = lept1.Eta();
    lept1_phi = lept1.Phi();

    deltaR_lept = lept0.DeltaR(lept1);

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

    if( cfg.selection()!="veryloose" && cfg.selection()!="presel" )
      if( gamma_pt/boss_mass< 40./150. ) continue;

    if( DATABLINDING && myTree.isData && boss_mass>500. ) continue;

    if( id==851 && doSyst ) { // systematic uncertainties

      // first scale
      float ref = myTree.LHEweight_original;

      float maxScaleDiff = 0.;

      TH1D* h1_pdf = new TH1D("pdf", "", 1000, -3000., 3000.);

      for( int i=0; i<myTree.nLHEweight; ++i ) {

        if( myTree.LHEweight_id[i]>=1000 && myTree.LHEweight_id[i]<1010 ) {
          float thisScaleDiff = fabs( (myTree.LHEweight_wgt[i]-ref)/ref );
          if( thisScaleDiff>maxScaleDiff) 
            maxScaleDiff = thisScaleDiff+1.;
        }

        if( myTree.LHEweight_id[i]>=2000 ) {
          if( myTree.LHEweight_wgt[i] > h1_pdf->GetXaxis()->GetXmax() || myTree.LHEweight_wgt[i] < h1_pdf->GetXaxis()->GetXmin() )
            std::cout << "WARNING!! PDF weight out of bounds: " << myTree.LHEweight_wgt[i] << std::endl;
          h1_pdf->Fill( myTree.LHEweight_wgt[i]/ref );
        }

      } // for lhe weights

      float meanPdfWgt = h1_pdf->GetMean();
      float rmsPdfWgt = h1_pdf->GetRMS();
      weight_pdf *= (1.+rmsPdfWgt/meanPdfWgt);

      weight_scale *= maxScaleDiff;

      delete h1_pdf;

    } // if correct sample

    outTree->Fill();
    
  } // for entries
    
  
  outTree->Write();

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
