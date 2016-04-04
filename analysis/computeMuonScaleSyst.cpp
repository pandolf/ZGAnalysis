#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TLorentzVector.h"

#include "../interface/rochcor2015.h"
#include "../interface/muresolution_run2.h"

#include "../interface/ZGDrawTools.h"
#include "../interface/ZGConfig.h"



int iComputedWithToys = 0;


void computeMuonSyst( const std::string& outdir, TFile* file, int mass, TGraph* gr_muSyst, TGraph* gr_egmSyst );
float computePtSyst( const TLorentzVector& lept, int pdgId, float sign );
float computePtSystFromToys( const TLorentzVector& lept, int pdgId );



int main( int argc, char* argv[] ) {


  if( argc<2 ) {
    std::cout << "USAGE: ./computeMuonScaleSyst [configFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  std::string configFileName(argv[1]);
  ZGConfig cfg(configFileName);

  std::string dir = cfg.getEventYieldDir() + "/muonScaleSyst";
  system(Form("mkdir -p %s", dir.c_str()));

  TFile* file = TFile::Open(Form("%s/trees.root", cfg.getEventYieldDir().c_str()));

  std::vector<int> masses;
  masses.push_back(350);
  masses.push_back(450);
  masses.push_back(600);
  masses.push_back(750);
  masses.push_back(900);
  masses.push_back(1500);
  masses.push_back(2000);

  TFile* outfile = TFile::Open(Form("%s/muonScaleSyst.root", dir.c_str()), "recreate" );

  TGraph* gr_muonScaleSyst = new TGraph(0);
  TGraph* gr_egmScaleSyst = new TGraph(0);

  for( unsigned i=0; i<masses.size(); ++i ) {
    std::cout << "-> Starting mass: " << masses[i] << std::endl;
    computeMuonSyst( dir, file, masses[i], gr_muonScaleSyst, gr_egmScaleSyst );
  }
  
  outfile->cd();
  gr_muonScaleSyst->SetName("gr_muScaleSyst");
  gr_egmScaleSyst->SetName("gr_egmScaleSyst");

  TF1* f1 = new TF1( "f1_muScaleSyst", "[0]/(1.+exp(-[1]*(x-[2])))", 200., 2100.);
  f1->SetParameter(2, 750.);
  f1->SetParameter(0, 0.05);
  gr_muonScaleSyst->Fit(f1, "RX");

  gr_muonScaleSyst->Write();
  gr_egmScaleSyst->Write();
  f1->Write();
  outfile->Close();

  std::cout << "-> Find your shit in: " << dir << std::endl;

  return 0;

}




void computeMuonSyst( const std::string& outdir, TFile* file, int mass, TGraph* gr_muSyst, TGraph* gr_egmSyst ) {


  iComputedWithToys = 0; // init for this mass point

  TTree* tree = (TTree*)file->Get(Form("XZg_Spin0_ZToLL_W_0p014_M_%d", mass) );



  TH1D* h1_bossMass_nominal = new TH1D("bossMass_nominal", "", 100, 0.5*(float)mass, 1.5*(float)mass );
  h1_bossMass_nominal->Sumw2();
  TH1D* h1_bossMass_muSystUp = new TH1D("bossMass_muSystUp", "", 100, 0.5*(float)mass, 1.5*(float)mass );
  h1_bossMass_muSystUp->Sumw2();
  TH1D* h1_bossMass_muSystDn = new TH1D("bossMass_muSystDn", "", 100, 0.5*(float)mass, 1.5*(float)mass );
  h1_bossMass_muSystDn->Sumw2();
  TH1D* h1_bossMass_egmSystUp = new TH1D("bossMass_egmSystUp", "", 100, 0.5*(float)mass, 1.5*(float)mass );
  h1_bossMass_egmSystUp->Sumw2();


  int leptType;
  tree->SetBranchAddress("leptType", &leptType);

  float boss_mass;
  tree->SetBranchAddress("boss_mass", &boss_mass);

  float lept0_pt;
  tree->SetBranchAddress("lept0_pt", &lept0_pt);
  float lept0_eta;
  tree->SetBranchAddress("lept0_eta", &lept0_eta);
  float lept0_phi;
  tree->SetBranchAddress("lept0_phi", &lept0_phi);
  float lept0_mass;
  tree->SetBranchAddress("lept0_mass", &lept0_mass);
  int lept0_pdgId;
  tree->SetBranchAddress("lept0_pdgId", &lept0_pdgId);

  float lept1_pt;
  tree->SetBranchAddress("lept1_pt", &lept1_pt);
  float lept1_eta;
  tree->SetBranchAddress("lept1_eta", &lept1_eta);
  float lept1_phi;
  tree->SetBranchAddress("lept1_phi", &lept1_phi);
  float lept1_mass;
  tree->SetBranchAddress("lept1_mass", &lept1_mass);
  int lept1_pdgId;
  tree->SetBranchAddress("lept1_pdgId", &lept1_pdgId);

  float gamma_pt;
  tree->SetBranchAddress("gamma_pt", &gamma_pt);
  float gamma_eta;
  tree->SetBranchAddress("gamma_eta", &gamma_eta);
  float gamma_phi;
  tree->SetBranchAddress("gamma_phi", &gamma_phi);
  float gamma_mass;
  tree->SetBranchAddress("gamma_mass", &gamma_mass);


  int nentries = tree->GetEntries();

  for( int iEntry=0; iEntry<nentries; ++iEntry ) {

    //std::cout << " Entry: " << iEntry << " / " << nentries << std::endl;
    if( iEntry % 1000 == 0 ) std::cout << " Entry: " << iEntry << " / " << nentries << std::endl;

    tree->GetEntry(iEntry);

    if( leptType!=13 ) continue;

    TLorentzVector lept0, lept1;
    lept0.SetPtEtaPhiM( lept0_pt, lept0_eta, lept0_phi, lept0_mass );
    lept1.SetPtEtaPhiM( lept1_pt, lept1_eta, lept1_phi, lept1_mass );

    float lept0_pt_systUp = computePtSyst( lept0, lept0_pdgId, +1. );
    float lept1_pt_systUp = computePtSyst( lept1, lept1_pdgId, +1. );
    float lept0_pt_systDn = computePtSyst( lept0, lept0_pdgId, -1. );
    float lept1_pt_systDn = computePtSyst( lept1, lept1_pdgId, -1. );

    if( lept0_pt_systUp<0. ) continue;
    if( lept1_pt_systUp<0. ) continue;
    if( lept0_pt_systDn<0. ) continue;
    if( lept1_pt_systDn<0. ) continue;

    TLorentzVector lept0_up, lept0_dn, lept1_up, lept1_dn;
    lept0_up.SetPtEtaPhiM( lept0_pt_systUp, lept0_eta, lept0_phi, lept0_mass );
    lept0_dn.SetPtEtaPhiM( lept0_pt_systDn, lept0_eta, lept0_phi, lept0_mass );
    lept1_up.SetPtEtaPhiM( lept1_pt_systUp, lept1_eta, lept1_phi, lept1_mass );
    lept1_dn.SetPtEtaPhiM( lept1_pt_systDn, lept1_eta, lept1_phi, lept1_mass );

    TLorentzVector gamma;
    gamma.SetPtEtaPhiM( gamma_pt, gamma_eta, gamma_phi, gamma_mass );

    TLorentzVector gamma_systUp;
    gamma_systUp.SetPtEtaPhiM( gamma_pt*1.01, gamma_eta, gamma_phi, gamma_mass );

    TLorentzVector boss_muSystUp = lept0_up + lept1_up + gamma;
    TLorentzVector boss_muSystDn = lept0_dn + lept1_dn + gamma;
    TLorentzVector boss_egmSystUp = lept0 + lept1 + gamma_systUp;

    h1_bossMass_nominal->Fill( boss_mass );
    h1_bossMass_muSystUp ->Fill( boss_muSystUp.M() );
    h1_bossMass_muSystDn ->Fill( boss_muSystDn.M() );
    h1_bossMass_egmSystUp->Fill( boss_egmSystUp.M() );

  } // for entries


  TFile* outfile = TFile::Open( Form("%s/muSyst_m%d.root", outdir.c_str(), mass), "recreate" );
  outfile->cd();

  h1_bossMass_muSystUp->SetLineColor(kGreen);
  h1_bossMass_muSystDn->SetLineColor(kRed);
  h1_bossMass_egmSystUp->SetLineColor(kRed);

  h1_bossMass_nominal->Write();
  h1_bossMass_muSystUp->Write();
  h1_bossMass_muSystDn->Write();
  h1_bossMass_egmSystUp->Write();

  outfile->Write();

  float muSystUp = fabs(h1_bossMass_muSystUp->GetMean()-h1_bossMass_nominal->GetMean())/h1_bossMass_nominal->GetMean();
  float muSystDn = fabs(h1_bossMass_muSystDn->GetMean()-h1_bossMass_nominal->GetMean())/h1_bossMass_nominal->GetMean();
  float egmSystUp = fabs(h1_bossMass_egmSystUp->GetMean()-h1_bossMass_nominal->GetMean())/h1_bossMass_nominal->GetMean();

  gr_muSyst ->SetPoint( gr_muSyst ->GetN(), mass, TMath::Max(muSystUp, muSystDn) );
  gr_egmSyst->SetPoint( gr_egmSyst->GetN(), mass, egmSystUp );

}





float computePtSyst( const TLorentzVector& lept, int pdgId, float sign ) {

  float returnPt = -1.;
  float pt = lept.Pt();

  if( pt < 200. ) {

    if( iComputedWithToys<2000 )
      returnPt = computePtSystFromToys(lept, pdgId);
    else
      returnPt = pt;
    iComputedWithToys++;

  } else {

    float scale = (fabs(lept.Eta())<1.479) ? 1.1 : 1.2;

    float ptTeV = pt/1000.;
    float oneOverPtTeV = 1./ptTeV;

    if( sign>0 )
      oneOverPtTeV *= scale;
    else
      oneOverPtTeV /= scale;

    returnPt = 1000./oneOverPtTeV;

  }

  return returnPt;

}



float computePtSystFromToys( const TLorentzVector& lept, int pdgId ) {

  //return 1.;

  int ntoys = 100;
  TH1D* h1_pt = new TH1D("pt", "", 100, 0.5*lept.Pt(), 1.5*lept.Pt() );

  for( int itoy=0; itoy<ntoys; ++itoy ) {

    TLorentzVector lept0_corr(lept);
    rochcor2015 *rmcor = new rochcor2015(itoy); // different seed
    float qter = 1.0;
    rmcor->momcor_mc(lept0_corr, pdgId/(abs(pdgId)), 0, qter);
    h1_pt->Fill( lept0_corr.Pt() );
    delete rmcor;

  }

  float returnSyst = lept.Pt()*(1.+h1_pt->GetRMS()/h1_pt->GetMean());
  delete h1_pt;

  return returnSyst;

}
