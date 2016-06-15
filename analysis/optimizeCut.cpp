#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"

#include "../interface/ZGDrawTools.h"



TH1D* optimizeCut( const std::string& outdir, TTree* tree_zg, TTree* tree_sig, float mass, float width, float crossSec, float eff_sel );



float lumi = 10.;




int main() {


  TFile* file = TFile::Open("~/CMSSW_7_6_3_ZG/src/ZGAnalysis/analysis/EventYields_presel_eth74X/trees.root");

  TTree* tree_zg = (TTree*)file->Get("zg");
  //TTree* tree_dy = (TTree*)file->Get("dy"); //ignore high-weight DY and multiply ZG by 1.1
  TTree* tree_375 = (TTree*)file->Get("XZg_Spin0_ZToLL_W_0p014_M_375");
  TTree* tree_750 = (TTree*)file->Get("XZg_Spin0_ZToLL_W_0p014_M_750");
  TTree* tree_1500 = (TTree*)file->Get("XZg_Spin0_ZToLL_W_0p014_M_1500");

  std::string outdir = "optCuts";
  system( Form("mkdir -p %s", outdir.c_str()) );


  TH1D* h1_sig_375 = optimizeCut( outdir, tree_zg, tree_375, 375., 0.03, 200., 0.5 );

  TFile* outfile = TFile::Open( Form("%s/file_opt.root", outdir.c_str()), "recreate");
  outfile->cd();

  h1_sig_375->Write();

  outfile->Close();
   
  return 0;

}


TH1D* optimizeCut( const std::string& outdir, TTree* tree_zg, TTree* tree_sig, float mass, float width, float crossSec,  float eff_sel ) {


  std::cout << "++ STARTING SIGNAL: " << mass << std::endl;

  float minCut = 0.;
  float maxCut = 0.7;
  float step = 0.02;
  int nSteps = (maxCut-minCut)/step;

  TH1D* h1_significance = new TH1D( Form("sig_%.0f", mass), "", nSteps, minCut, maxCut );

  TH1D* h1_ptOm_zg  = new TH1D("ptOm_zg" , "", nSteps, minCut, maxCut );
  TH1D* h1_ptOm_sig = new TH1D("ptOm_sig", "", nSteps, minCut, maxCut );

  h1_ptOm_zg ->Sumw2();
  h1_ptOm_sig->Sumw2();

  tree_zg ->Project( "ptOm_zg" , "gamma_pt/boss_mass", "1.1*weight");
  tree_sig->Project( "ptOm_sig", "gamma_pt/boss_mass",     "weight");


  TH1D* h1_ptOm_massCut_zg  = new TH1D("ptOm_massCut_zg" , "", nSteps, minCut, maxCut );
  TH1D* h1_ptOm_massCut_sig = new TH1D("ptOm_massCut_sig", "", nSteps, minCut, maxCut );

  h1_ptOm_massCut_zg ->Sumw2();
  h1_ptOm_massCut_sig->Sumw2();


  float massMin = mass*(1.-width);
  float massMax = mass*(1.+width);

  std::cout << "  mass cuts: " << massMin << "-" << massMax << std::endl;

  tree_zg ->Project( "ptOm_massCut_zg" , "gamma_pt/boss_mass", Form("1.1*weight*(boss_mass>%f && boss_mass<%f)", massMin, massMax) );
  tree_sig->Project( "ptOm_massCut_sig", "gamma_pt/boss_mass", Form(    "weight*(boss_mass>%f && boss_mass<%f)", massMin, massMax) );


  float effMassCutSig = h1_ptOm_massCut_sig->Integral()/h1_ptOm_sig->Integral();

  std::cout << "  efficiency of mass cut: " << effMassCutSig << std::endl;
  std::cout << std::endl;


  //for( float thresh = minCut; thresh += step; cut<=maxCut ) { 
  for( int iStep=0; iStep<nSteps; ++iStep ) {

    float thresh = minCut + (float)iStep*step;

    std::cout << "step " << iStep << "::   cut: " << thresh << std::endl;

    int binMin = h1_ptOm_massCut_zg->FindBin( thresh );
    int binMax = h1_ptOm_massCut_zg->FindBin( maxCut );

    float thisBG = h1_ptOm_massCut_zg->Integral( binMin, binMax );

    float thisCutEffSig = h1_ptOm_massCut_sig->Integral( binMin, binMax )/h1_ptOm_massCut_sig->Integral();

    float thisSig = crossSec * eff_sel * thisCutEffSig * effMassCutSig * 2.*0.0033;

    thisBG  *= lumi;
    thisSig *= lumi;

    std::cout << " bg: " << thisBG << std::endl;
    std::cout << " sig: " << thisSig << std::endl;

    if( thisBG <=0. ) continue;

    h1_significance->SetBinContent( binMin, thisSig/sqrt(thisBG) );

  }

  return h1_significance;

}
