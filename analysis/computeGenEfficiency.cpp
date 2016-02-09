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


#include "interface/ZGSample.h"
#include "interface/ZGCommonTools.h"



#define zg_cxx
#include "interface/zg.h"




int main( int argc, char* argv[] ) {




  std::string filename = "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/mmasciov/MT2production/74X/Spring15/PostProcessed/20Dec2015_forGunther/ZGTo2LG_post.root";
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



  TFile* outFile = TFile::Open("provaGenEfficiency.root", "recreate");
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

  TH1D* h1_eff_all_num = new TH1D( "eff_all_num", "", nBins, bins );
  h1_eff_all_num->Sumw2();
  TH1D* h1_eff_noHLT_num = new TH1D( "eff_noHLT_num", "", nBins, bins );
  h1_eff_noHLT_num->Sumw2();
  TH1D* h1_eff_noIso_num = new TH1D( "eff_noIso_num", "", nBins, bins );
  h1_eff_noIso_num->Sumw2();

  TH1D* h1_massBias = new TH1D( "massBias", "", nBins, bins );
  h1_massBias->Sumw2();
  TH1D* h1_massReso = new TH1D( "massReso", "", nBins, bins );
  h1_massReso->Sumw2();

  std::vector<TH1D*> vh1_massReso;
  for( int i=0; i<nBins; ++i ) {
    TH1D* h1_tmp = new TH1D( Form("reso_%d", i), "", 40, -0.5, 0.5);
    h1_tmp->Sumw2();
    vh1_massReso.push_back( h1_tmp );
  }


 
  int nentries = tree->GetEntries();
  
  for( int iEntry=0; iEntry<nentries; ++iEntry ) {
    
    if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;
    
    myTree.GetEntry(iEntry);

    float weight = (myTree.isData) ? 1. : myTree.evt_scale1fb;
    // pu reweighting:
    if( !myTree.isData ) {
      float puWeight = ZGCommonTools::getPUweight( myTree.nVert, h1_nVert_data, h1_nVert_mc );
      weight *= puWeight;
    }




    // first find leptons
    if( myTree.ngenLep!=2 ) continue;
    TLorentzVector genLep0, genLep1;
    genLep0.SetPtEtaPhiM( myTree.genLep_pt[0], myTree.genLep_eta[0], myTree.genLep_phi[0], myTree.genLep_mass[0] );
    genLep1.SetPtEtaPhiM( myTree.genLep_pt[1], myTree.genLep_eta[1], myTree.genLep_phi[1], myTree.genLep_mass[1] );


    float maxPt = TMath::Max( genLep0.Pt(), genLep1.Pt() );
    float minPt = TMath::Min( genLep0.Pt(), genLep1.Pt() );
    if( maxPt<25. ) continue;
    if( minPt<20. ) continue;

    if( fabs(genLep0.Eta()) > 2.4 ) continue;
    if( fabs(genLep1.Eta()) > 2.4 ) continue;

    TLorentzVector genZ = genLep0 + genLep1;
    if( genZ.M()<50. ) continue;


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


    // and now reco!

    if( myTree.nVert==0 ) continue;
    
    
    if( myTree.nlep!=2 ) continue; // two leptons
    if( myTree.lep_pdgId[0] != -myTree.lep_pdgId[1] ) continue; // same flavour, opposite sign

    int leptType = abs(myTree.lep_pdgId[0]);
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

    TLorentzVector zBoson = lept0+lept1;
    if( zBoson.M()<50. ) continue;


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
      float deltaR_thresh = 0.4;
      if( tmp_photon.DeltaR(lept0)<deltaR_thresh || tmp_photon.DeltaR(lept1)<deltaR_thresh ) continue;

      photon.SetPtEtaPhiM( myTree.gamma_pt[iPhot], myTree.gamma_eta[iPhot], myTree.gamma_phi[iPhot], myTree.gamma_mass[iPhot] );
      foundPhoton = true;

    }

    if( !foundPhoton ) continue;


    TLorentzVector boss = zBoson + photon;
    float recoMass = boss.M();
    if( recoMass<200. ) continue;

    h1_eff_noHLT_num->Fill( genMass, weight );

    if( !( myTree.HLT_DoubleEl || myTree.HLT_DoubleMu || myTree.HLT_Photon165_HE10)) continue;

    h1_eff_noIso_num->Fill( genMass, weight );

    if( myTree.gamma_chHadIso[0]>2.5 ) continue;

    h1_eff_all_num->Fill( genMass, weight );

    int iBin = h1_eff_all_num->FindBin( genMass ) - 1;

    if( iBin>=0 && iBin<vh1_massReso.size() )
      vh1_massReso[iBin]->Fill( (recoMass-genMass)/genMass, weight );

    
  } // for entries
    

  outFile->cd();

  h1_eff_denom->Write();
  h1_eff_all_num->Write();
  h1_eff_noHLT_num->Write();
  h1_eff_noIso_num->Write();

  TEfficiency* eff_all = new TEfficiency( *h1_eff_all_num, *h1_eff_denom);
  eff_all->SetName( "eff_all" );
  eff_all->Write();

  TEfficiency* eff_noHLT = new TEfficiency( *h1_eff_noHLT_num, *h1_eff_denom);
  eff_noHLT->SetName( "eff_noHLT" );
  eff_noHLT->Write();

  TEfficiency* eff_noIso = new TEfficiency( *h1_eff_noIso_num, *h1_eff_denom);
  eff_noIso->SetName( "eff_noIso" );
  eff_noIso->Write();

  
  for( unsigned i=0; i<vh1_massReso.size(); ++i ) {
    h1_massBias->SetBinContent( i+1, vh1_massReso[i]->GetMean() );
    h1_massReso->SetBinContent( i+1, vh1_massReso[i]->GetRMS() );
    h1_massBias->SetBinError( i+1, vh1_massReso[i]->GetMeanError() );
    h1_massReso->SetBinError( i+1, vh1_massReso[i]->GetRMSError() );
    vh1_massReso[i]->Write();
  }

  h1_massReso->Write();
  h1_massBias->Write();
  
  outFile->Close();

  return 0;

}


