#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>


#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TEfficiency.h"
#include "TRandom3.h"


#include "interface/ZGSample.h"
#include "interface/ZGCommonTools.h"

// muon rochester corrections:
#include "../interface/rochcor2015.h"
#include "../interface/muresolution_run2.h"


#define zg_cxx
#include "interface/zg.h"


TRandom3 myRandom_(13);

void smearEmEnergy( TLorentzVector& p );


bool use76 = false;


int main( int argc, char* argv[] ) {


  std::string filename;
  if( use76 ) 
    filename = "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/mschoene/MT2production/76X/Spring15/PostProcessed/12Feb2016_ZGammaMC/ZGTo2LG_post.root";
  else
    filename = "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/mmasciov/MT2production/74X/Spring15/PostProcessed/20Dec2015_forGunther/ZGTo2LG_post.root";
  TFile* file = TFile::Open(filename.c_str());
  TTree* tree = (TTree*)file->Get("mt2");

  std::cout << "-> Opened file: " << filename << std::endl;



  TFile* puFile_data = TFile::Open("puData.root");
  TH1D* h1_nVert_data = (TH1D*)puFile_data->Get("nVert");

  TFile* puFile_mc = TFile::Open("puMC.root");
  TH1D* h1_nVert_mc = (TH1D*)puFile_mc->Get("nVert");


  ZGTree myTree;
  myTree.loadGenStuff = true;
  myTree.Init(tree);



  TFile* outFile;
  if( use76 )
    outFile = TFile::Open("genEfficiency76.root", "recreate");
  else
    outFile = TFile::Open("genEfficiency.root", "recreate");
  outFile->cd();

  int nBins = 8;
  Double_t bins[nBins+1];
  bins[0] = 300.;
  bins[1] = 350.;
  bins[2] = 400.;
  bins[3] = 450.;
  bins[4] = 500.;
  bins[5] = 600.;
  bins[6] = 700.;
  bins[7] = 800.;
  bins[8] = 950.;
  
  TH1D* h1_eff_denom = new TH1D( "eff_denom", "", nBins, bins );
  h1_eff_denom->Sumw2();
  TH1D* h1_eff_denom_ee = new TH1D( "eff_denom_ee", "", nBins, bins );
  h1_eff_denom_ee->Sumw2();
  TH1D* h1_eff_denom_mm = new TH1D( "eff_denom_mm", "", nBins, bins );
  h1_eff_denom_mm->Sumw2();

  TH1D* h1_eff_all_num = new TH1D( "eff_all_num", "", nBins, bins );
  h1_eff_all_num->Sumw2();
  TH1D* h1_eff_ee_num = new TH1D( "eff_ee_num", "", nBins, bins );
  h1_eff_ee_num->Sumw2();
  TH1D* h1_eff_mm_num = new TH1D( "eff_mm_num", "", nBins, bins );
  h1_eff_mm_num->Sumw2();
  TH1D* h1_eff_noHLT_num = new TH1D( "eff_noHLT_num", "", nBins, bins );
  h1_eff_noHLT_num->Sumw2();
  TH1D* h1_eff_noIso_num = new TH1D( "eff_noIso_num", "", nBins, bins );
  h1_eff_noIso_num->Sumw2();

  TH1D* h1_massBias = new TH1D( "massBias", "", nBins, bins );
  h1_massBias->Sumw2();
  TH1D* h1_massReso = new TH1D( "massReso", "", nBins, bins );
  h1_massReso->Sumw2();

  TH1D* h1_massBias_ee = new TH1D( "massBias_ee", "", nBins, bins );
  h1_massBias_ee->Sumw2();
  TH1D* h1_massReso_ee = new TH1D( "massReso_ee", "", nBins, bins );
  h1_massReso_ee->Sumw2();

  TH1D* h1_massBias_mm = new TH1D( "massBias_mm", "", nBins, bins );
  h1_massBias_mm->Sumw2();
  TH1D* h1_massReso_mm = new TH1D( "massReso_mm", "", nBins, bins );
  h1_massReso_mm->Sumw2();

  std::vector<TH1D*> vh1_massReso, vh1_massReso_ee, vh1_massReso_mm;

  for( int i=0; i<nBins; ++i ) {

    TH1D* h1_tmp = new TH1D( Form("reso_%d", i), "", 40, -0.2, 0.2);
    h1_tmp->Sumw2();
    vh1_massReso.push_back( h1_tmp );

    TH1D* h1_tmp_ee = new TH1D( Form("reso_ee_%d", i), "", 40, -0.2, 0.2);
    h1_tmp_ee->Sumw2();
    vh1_massReso_ee.push_back( h1_tmp_ee );

    TH1D* h1_tmp_mm = new TH1D( Form("reso_mm_%d", i), "", 40, -0.2, 0.2);
    h1_tmp_mm->Sumw2();
    vh1_massReso_mm.push_back( h1_tmp_mm );

  }

  rochcor2015 *rmcor = new rochcor2015();


  TTree* outtree = new TTree("genTree", "");

  int leptType;
  outtree->Branch( "leptType", &leptType, "leptType/I" );

  float gammaReco_pt;
  outtree->Branch( "gammaReco_pt", &gammaReco_pt, "gammaReco_pt/F" );
  float gammaReco_eta;
  outtree->Branch( "gammaReco_eta", &gammaReco_eta, "gammaReco_eta/F" );
  float gammaReco_phi;
  outtree->Branch( "gammaReco_phi", &gammaReco_phi, "gammaReco_phi/F" );
  float gammaReco_mass;
  outtree->Branch( "gammaReco_mass", &gammaReco_mass, "gammaReco_mass/F" );

  float lept0Reco_pt;
  outtree->Branch( "lept0Reco_pt", &lept0Reco_pt, "lept0Reco_pt/F" );
  float lept0Reco_eta;
  outtree->Branch( "lept0Reco_eta", &lept0Reco_eta, "lept0Reco_eta/F" );
  float lept0Reco_phi;
  outtree->Branch( "lept0Reco_phi", &lept0Reco_phi, "lept0Reco_phi/F" );
  float lept0Reco_mass;
  outtree->Branch( "lept0Reco_mass", &lept0Reco_mass, "lept0Reco_mass/F" );

  float lept1Reco_pt;
  outtree->Branch( "lept1Reco_pt", &lept1Reco_pt, "lept1Reco_pt/F" );
  float lept1Reco_eta;
  outtree->Branch( "lept1Reco_eta", &lept1Reco_eta, "lept1Reco_eta/F" );
  float lept1Reco_phi;
  outtree->Branch( "lept1Reco_phi", &lept1Reco_phi, "lept1Reco_phi/F" );
  float lept1Reco_mass;
  outtree->Branch( "lept1Reco_mass", &lept1Reco_mass, "lept1Reco_mass/F" );

  float zReco_pt;
  outtree->Branch( "zReco_pt", &zReco_pt, "zReco_pt/F" );
  float zReco_eta;
  outtree->Branch( "zReco_eta", &zReco_eta, "zReco_eta/F" );
  float zReco_phi;
  outtree->Branch( "zReco_phi", &zReco_phi, "zReco_phi/F" );
  float zReco_mass;
  outtree->Branch( "zReco_mass", &zReco_mass, "zReco_mass/F" );

  float bossReco_pt;
  outtree->Branch( "bossReco_pt", &bossReco_pt, "bossReco_pt/F" );
  float bossReco_eta;
  outtree->Branch( "bossReco_eta", &bossReco_eta, "bossReco_eta/F" );
  float bossReco_phi;
  outtree->Branch( "bossReco_phi", &bossReco_phi, "bossReco_phi/F" );
  float bossReco_mass;
  outtree->Branch( "bossReco_mass", &bossReco_mass, "bossReco_mass/F" );

  float gammaGen_pt;
  outtree->Branch( "gammaGen_pt", &gammaGen_pt, "gammaGen_pt/F" );
  float gammaGen_eta;
  outtree->Branch( "gammaGen_eta", &gammaGen_eta, "gammaGen_eta/F" );
  float gammaGen_phi;
  outtree->Branch( "gammaGen_phi", &gammaGen_phi, "gammaGen_phi/F" );
  float gammaGen_mass;
  outtree->Branch( "gammaGen_mass", &gammaGen_mass, "gammaGen_mass/F" );

  float lept0Gen_pt;
  outtree->Branch( "lept0Gen_pt", &lept0Gen_pt, "lept0Gen_pt/F" );
  float lept0Gen_eta;
  outtree->Branch( "lept0Gen_eta", &lept0Gen_eta, "lept0Gen_eta/F" );
  float lept0Gen_phi;
  outtree->Branch( "lept0Gen_phi", &lept0Gen_phi, "lept0Gen_phi/F" );
  float lept0Gen_mass;
  outtree->Branch( "lept0Gen_mass", &lept0Gen_mass, "lept0Gen_mass/F" );

  float lept1Gen_pt;
  outtree->Branch( "lept1Gen_pt", &lept1Gen_pt, "lept1Gen_pt/F" );
  float lept1Gen_eta;
  outtree->Branch( "lept1Gen_eta", &lept1Gen_eta, "lept1Gen_eta/F" );
  float lept1Gen_phi;
  outtree->Branch( "lept1Gen_phi", &lept1Gen_phi, "lept1Gen_phi/F" );
  float lept1Gen_mass;
  outtree->Branch( "lept1Gen_mass", &lept1Gen_mass, "lept1Gen_mass/F" );

  float zGen_pt;
  outtree->Branch( "zGen_pt", &zGen_pt, "zGen_pt/F" );
  float zGen_eta;
  outtree->Branch( "zGen_eta", &zGen_eta, "zGen_eta/F" );
  float zGen_phi;
  outtree->Branch( "zGen_phi", &zGen_phi, "zGen_phi/F" );
  float zGen_mass;
  outtree->Branch( "zGen_mass", &zGen_mass, "zGen_mass/F" );

  float bossGen_pt;
  outtree->Branch( "bossGen_pt", &bossGen_pt, "bossGen_pt/F" );
  float bossGen_eta;
  outtree->Branch( "bossGen_eta", &bossGen_eta, "bossGen_eta/F" );
  float bossGen_phi;
  outtree->Branch( "bossGen_phi", &bossGen_phi, "bossGen_phi/F" );
  float bossGen_mass;
  outtree->Branch( "bossGen_mass", &bossGen_mass, "bossGen_mass/F" );


 
  int nentries = tree->GetEntries();
  
  for( int iEntry=0; iEntry<nentries; ++iEntry ) {
    
    if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;
    
    myTree.GetEntry(iEntry);

    float weight = (myTree.isData) ? 1. : myTree.evt_scale1fb;
    // pu reweighting:
    if( !myTree.isData ) {
      //weight *= myTree.puWeight;
    }


    // first find leptons
    //if( myTree.ngenLep!=2 ) continue;
    //TLorentzVector genLep0, genLep1;
    //genLep0.SetPtEtaPhiM( myTree.genLep_pt[0], myTree.genLep_eta[0], myTree.genLep_phi[0], myTree.genLep_mass[0] );
    //genLep1.SetPtEtaPhiM( myTree.genLep_pt[1], myTree.genLep_eta[1], myTree.genLep_phi[1], myTree.genLep_mass[1] );
    
    std::vector<TLorentzVector> genLeptons;
    int genLeptType = 0;
    for( int iGen=0; iGen<myTree.ngenPart && genLeptons.size()<2; ++iGen ) {

      if( myTree.genPart_pt[iGen]<1. ) continue;
      if( myTree.genPart_status[iGen]!=1 ) continue;
      if( abs(myTree.genPart_pdgId[iGen])!=11 && abs(myTree.genPart_pdgId[iGen])!=13 ) continue;
      if( myTree.genPart_motherId[iGen]!=myTree.genPart_pdgId[iGen] ) continue;
      TLorentzVector tmpLep;
      tmpLep.SetPtEtaPhiM( myTree.genPart_pt[iGen], myTree.genPart_eta[iGen], myTree.genPart_phi[iGen], myTree.genPart_mass[iGen] );
      genLeptons.push_back(tmpLep);
      genLeptType = abs(myTree.genPart_pdgId[iGen]);

    }

    if( genLeptType!=11 && genLeptType!=13 ) continue;
    if( genLeptons.size()<2 ) continue;

    TLorentzVector genLep0, genLep1;
    genLep0 = genLeptons[0];
    genLep1 = genLeptons[1];

    float maxPt = TMath::Max( genLep0.Pt(), genLep1.Pt() );
    float minPt = TMath::Min( genLep0.Pt(), genLep1.Pt() );
    if( maxPt<25. ) continue;
    if( minPt<20. ) continue;

    if( fabs(genLep0.Eta()) > 2.4 ) continue;
    if( fabs(genLep1.Eta()) > 2.4 ) continue;

    TLorentzVector genZ = genLep0 + genLep1;
    if( genZ.M()<50. ) continue;
    //if( genZ.M()<50. || genZ.M()>130. ) continue;


    TLorentzVector genPhoton;
    bool foundGenPhoton = false;

    for( int iGen=0; iGen<myTree.ngenPart && !foundGenPhoton; ++iGen ) {

      if( myTree.genPart_pdgId[iGen]!=22 ) continue;
      if( myTree.genPart_status[iGen]!=1 ) continue;
      if( myTree.genPart_pt[iGen]<40. ) continue;
      if( fabs(myTree.genPart_eta[iGen])>2.5 ) continue;

      TLorentzVector photon_temp;
      photon_temp.SetPtEtaPhiM( myTree.genPart_pt[iGen], myTree.genPart_eta[iGen], myTree.genPart_phi[iGen], myTree.genPart_mass[iGen] );

      float deltaRmin_part = 999.;
      // look for closest parton
      for( int jGen=0; jGen<myTree.ngenPart && iGen!=jGen; ++jGen ) {

        if( myTree.genPart_pt[jGen]<1. ) continue;
        if( myTree.genPart_status[jGen]<=21 ) continue;
        bool isParton = ( myTree.genPart_pdgId[jGen]==21 || abs(myTree.genPart_pdgId[jGen])<7 );
        if( !isParton ) continue;

        TLorentzVector thisparton;
        thisparton.SetPtEtaPhiM( myTree.genPart_pt[jGen], myTree.genPart_eta[jGen], myTree.genPart_phi[jGen], myTree.genPart_mass[jGen] );

        float thisDeltaR = thisparton.DeltaR(photon_temp);
        if( thisDeltaR<deltaRmin_part ) {
          deltaRmin_part = thisDeltaR;
        }
    
      }

      // far away from partons
      if( deltaRmin_part<0.4 ) continue;


      // far away from leptons
      if( photon_temp.DeltaR( genLep0 ) > 0.4 && photon_temp.DeltaR( genLep1 ) > 0.4 ) {
        genPhoton.SetPtEtaPhiM( myTree.genPart_pt[iGen], myTree.genPart_eta[iGen], myTree.genPart_phi[iGen], myTree.genPart_mass[iGen] );
        foundGenPhoton = true;
      }

    }


    if( !foundGenPhoton ) continue;
    if( genPhoton.Pt()<40. ) continue;
    if( fabs(genPhoton.Eta())>2.5 ) continue;
    if( fabs(genPhoton.Eta())>1.44 && fabs(genPhoton.Eta())<1.57 ) continue;


    TLorentzVector genBoss = genZ + genPhoton;
    float genMass = genBoss.M();
    if( genMass<200. ) continue;

    h1_eff_denom->Fill( genMass, weight );
    if( genLeptType==11 )
      h1_eff_denom_ee->Fill( genMass, weight );
    else
      h1_eff_denom_mm->Fill( genMass, weight );


    // and now reco!

    if( myTree.nVert==0 ) continue;
    
    
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

    if( leptType==11 ) { //electrons
      if( myTree.lep_tightId[0]==0 || myTree.lep_tightId[1]==0 ) continue; // loose electron ID
    } else { // muons
      float qter = 1.0;
      rmcor->momcor_mc(lept0, myTree.lep_pdgId[0]/(abs(myTree.lep_pdgId[0])), 0, qter);
      rmcor->momcor_mc(lept1, myTree.lep_pdgId[1]/(abs(myTree.lep_pdgId[1])), 0, qter);
      //if( myTree.lep_tightId[0]==0 && myTree.lep_tightId[1]==0 ) continue; // tight muon ID on one leg
    }


    TLorentzVector zBoson = lept0+lept1;
    if( zBoson.M()<50. ) continue;
    //if( zBoson.M()<50. || zBoson.M()>130. ) continue;


    if( myTree.ngamma==0 ) continue; // photon

    bool foundPhoton = false;
    TLorentzVector photon;

    for( int iPhot=0; iPhot<myTree.ngamma && !foundPhoton; ++iPhot ) {

      TLorentzVector tmp_photon;
      tmp_photon.SetPtEtaPhiM( myTree.gamma_pt[iPhot], myTree.gamma_eta[iPhot], myTree.gamma_phi[iPhot], myTree.gamma_mass[iPhot] );
 
      if( tmp_photon.Pt()<40. ) continue;
      if( fabs(tmp_photon.Eta())>1.44 && fabs(tmp_photon.Eta())<1.57 ) continue;
      if( fabs(tmp_photon.Eta())>2.5 ) continue;
      if( myTree.gamma_idCutBased[iPhot]==0 ) continue;
      //if( fabs(myTree.gamma_eta[iPhot])<1.44 ) {
      //  if( myTree.gamma_sigmaIetaIeta[iPhot]>0.0102 ) continue;
      //} else {
      //  if( myTree.gamma_sigmaIetaIeta[iPhot]>0.0274 ) continue;
      //}
      float deltaR_thresh = 0.4;
      if( tmp_photon.DeltaR(lept0)<deltaR_thresh || tmp_photon.DeltaR(lept1)<deltaR_thresh ) continue;

      photon.SetPtEtaPhiM( myTree.gamma_pt[iPhot], myTree.gamma_eta[iPhot], myTree.gamma_phi[iPhot], myTree.gamma_mass[iPhot] );
      foundPhoton = true;

    }

    if( !foundPhoton ) continue;

    smearEmEnergy( photon );

    TLorentzVector boss = zBoson + photon;
    float recoMass = boss.M();
    if( recoMass<200. ) continue;

    h1_eff_noHLT_num->Fill( genMass, weight );

    if( !( myTree.HLT_DoubleEl || myTree.HLT_DoubleMu || myTree.HLT_DoubleEl33 || myTree.HLT_SingleMu ) ) continue;

    h1_eff_noIso_num->Fill( genMass, weight );

    if( myTree.gamma_chHadIso[0]>2.5 ) continue;

    h1_eff_all_num->Fill( genMass, weight );
    if( leptType==11 ) h1_eff_ee_num->Fill( genMass, weight );
    else               h1_eff_mm_num->Fill( genMass, weight );

    int iBin = h1_eff_all_num->FindBin( genMass ) - 1;

    if( iBin>=0 && iBin<vh1_massReso.size() ) {
      vh1_massReso[iBin]->Fill( (recoMass-genMass)/genMass, weight );
      if( leptType==11 ) {
        vh1_massReso_ee[iBin]->Fill( (recoMass-genMass)/genMass, weight );
      } else {
        vh1_massReso_mm[iBin]->Fill( (recoMass-genMass)/genMass, weight );
      }
    }

    // set tree branches and fill tree
    gammaReco_pt   = photon.Pt();
    gammaReco_eta  = photon.Eta();
    gammaReco_phi  = photon.Phi();
    gammaReco_mass = photon.M();

    lept0Reco_pt   = lept0.Pt();
    lept0Reco_eta  = lept0.Eta();
    lept0Reco_phi  = lept0.Phi();
    lept0Reco_mass = lept0.M();
    
    lept1Reco_pt   = lept1.Pt();
    lept1Reco_eta  = lept1.Eta();
    lept1Reco_phi  = lept1.Phi();
    lept1Reco_mass = lept1.M();
    
    zReco_pt   = zBoson.Pt();
    zReco_eta  = zBoson.Eta();
    zReco_phi  = zBoson.Phi();
    zReco_mass = zBoson.M();
    
    bossReco_pt   = boss.Pt();
    bossReco_eta  = boss.Eta();
    bossReco_phi  = boss.Phi();
    bossReco_mass = boss.M();
    
    gammaGen_pt   = genPhoton.Pt();
    gammaGen_eta  = genPhoton.Eta();
    gammaGen_phi  = genPhoton.Phi();
    gammaGen_mass = genPhoton.M();

    lept0Gen_pt   = genLep0.Pt();
    lept0Gen_eta  = genLep0.Eta();
    lept0Gen_phi  = genLep0.Phi();
    lept0Gen_mass = genLep0.M();
    
    lept1Gen_pt   = genLep1.Pt();
    lept1Gen_eta  = genLep1.Eta();
    lept1Gen_phi  = genLep1.Phi();
    lept1Gen_mass = genLep1.M();
    
    zGen_pt   = genZ.Pt();
    zGen_eta  = genZ.Eta();
    zGen_phi  = genZ.Phi();
    zGen_mass = genZ.M();
    
    bossGen_pt   = genBoss.Pt();
    bossGen_eta  = genBoss.Eta();
    bossGen_phi  = genBoss.Phi();
    bossGen_mass = genBoss.M();

    outtree->Fill();
    
  } // for entries
    

  outFile->cd();

  h1_eff_denom->Write();
  h1_eff_denom_ee->Write();
  h1_eff_denom_mm->Write();
  h1_eff_all_num->Write();
  h1_eff_ee_num->Write();
  h1_eff_mm_num->Write();
  h1_eff_noHLT_num->Write();
  h1_eff_noIso_num->Write();

  TEfficiency* eff_all = new TEfficiency( *h1_eff_all_num, *h1_eff_denom);
  eff_all->SetName( "eff_all" );
  eff_all->Write();

  TEfficiency* eff_ee = new TEfficiency( *h1_eff_ee_num, *h1_eff_denom_ee);
  eff_ee->SetName( "eff_ee" );
  eff_ee->Write();

  TEfficiency* eff_mm = new TEfficiency( *h1_eff_mm_num, *h1_eff_denom_mm);
  eff_mm->SetName( "eff_mm" );
  eff_mm->Write();

  TEfficiency* eff_noHLT = new TEfficiency( *h1_eff_noHLT_num, *h1_eff_denom);
  eff_noHLT->SetName( "eff_noHLT" );
  eff_noHLT->Write();

  TEfficiency* eff_noIso = new TEfficiency( *h1_eff_noIso_num, *h1_eff_denom);
  eff_noIso->SetName( "eff_noIso" );
  eff_noIso->Write();

  
  for( unsigned i=0; i<vh1_massReso.size(); ++i ) {

    h1_massBias->SetBinContent( i+1, vh1_massReso[i]->GetMean() );
    h1_massReso->SetBinContent( i+1, vh1_massReso[i]->GetRMS() );
    h1_massBias->SetBinError  ( i+1, vh1_massReso[i]->GetMeanError() );
    h1_massReso->SetBinError  ( i+1, vh1_massReso[i]->GetRMSError() );

    h1_massBias_ee->SetBinContent( i+1, vh1_massReso_ee[i]->GetMean() );
    h1_massReso_ee->SetBinContent( i+1, vh1_massReso_ee[i]->GetRMS() );
    h1_massBias_ee->SetBinError  ( i+1, vh1_massReso_ee[i]->GetMeanError() );
    h1_massReso_ee->SetBinError  ( i+1, vh1_massReso_ee[i]->GetRMSError() );

    h1_massBias_mm->SetBinContent( i+1, vh1_massReso_mm[i]->GetMean() );
    h1_massReso_mm->SetBinContent( i+1, vh1_massReso_mm[i]->GetRMS() );
    h1_massBias_mm->SetBinError  ( i+1, vh1_massReso_mm[i]->GetMeanError() );
    h1_massReso_mm->SetBinError  ( i+1, vh1_massReso_mm[i]->GetRMSError() );

    vh1_massReso[i]->Write();
    vh1_massReso_ee[i]->Write();
    vh1_massReso_mm[i]->Write();

  }

  h1_massReso->Write();
  h1_massBias->Write();
  h1_massReso_ee->Write();
  h1_massBias_ee->Write();
  h1_massReso_mm->Write();
  h1_massBias_mm->Write();

  outtree->Write();
  
  outFile->Close();

  return 0;

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


