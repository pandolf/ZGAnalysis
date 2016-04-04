#include "TH2D.h"
#include "TH1D.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"


#include "../interface/ZGDrawTools.h"




int main() {

  ZGDrawTools::setStyle();

  TFile* file = TFile::Open("EventYields_v0_eth74X_miniIso999/trees.root");
  //TFile* file = TFile::Open("ZG_miniIsoLt999.root");
  TTree* tree = (TTree*)file->Get("zg");
  //TTree* tree = (TTree*)file->Get("mt2");

  int nBins(14);
  int nBinsPlusOne(nBins+1);
  float xMin = 100.;
  float xMax = 1100.;

  Double_t bins[nBinsPlusOne];
  bins [0] = xMin;
  bins [1] = 120.; 
  bins [2] = 140.; 
  bins [3] = 160.; 
  bins [4] = 180.; 
  bins [5] = 200.; 
  bins [6] = 240.; 
  bins [7] = 280.; 
  bins [8] = 320.; 
  bins [9] = 360.; 
  bins[10] = 400.; 
  bins[11] = 500.; 
  bins[12] = 600.; 
  bins[13] = 800.; 
  bins[14] = xMax;

  TH1D* h1_num_ee = new TH1D("num_ee", "", nBins, bins);
  TH1D* h1_denom_ee = new TH1D("denom_ee", "", nBins, bins);

  //tree->Project("denom_ee", "lep_pt[0]", "abs(lep_pdgId[0])==11");
  //tree->Project("num_ee"  , "lep_pt[0]", "abs(lep_pdgId[0])==11 && lep_miniRelIso[0]<0.1");
  tree->Project("denom_ee", "boss_mass", "leptType==11");
  tree->Project("num_ee"  , "boss_mass", "leptType==11 && lept0_miniRelIso<0.1 && lept1_miniRelIso<0.1");
  
  TEfficiency* eff_ee = new TEfficiency( *h1_num_ee, *h1_denom_ee );

  TH1D* h1_num_mm = new TH1D("num_mm", "", nBins, bins);
  TH1D* h1_denom_mm = new TH1D("denom_mm", "", nBins, bins);

  //tree->Project("denom_mm", "lep_pt[0]", "abs(lep_pdgId[0])==13");
  //tree->Project("num_mm"  , "lep_pt[0]", "abs(lep_pdgId[0])==13 && lep_miniRelIso[0]<0.2");
  tree->Project("denom_mm", "boss_mass", "leptType==13");
  tree->Project("num_mm"  , "boss_mass", "leptType==13 && lept0_miniRelIso<0.2 && lept1_miniRelIso<0.2");
  
  TEfficiency* eff_mm = new TEfficiency( *h1_num_mm, *h1_denom_mm );


  TCanvas* c1 = new TCanvas("c1", "", 600, 600 );
  c1->cd();

  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0.7, 1.001 );
  h2_axes->SetXTitle( "M(Z#gamma) [GeV]");
  //h2_axes->SetXTitle( "Lepton p_{T} [GeV]"); 
  h2_axes->SetYTitle( "Mini Iso Cut Efficiency" );
  h2_axes->Draw();



  eff_ee->SetMarkerStyle(20);
  eff_ee->SetMarkerSize(1.5);
  eff_ee->SetMarkerColor(kRed);
  eff_ee->SetLineColor(kRed);

  eff_mm->SetMarkerStyle(20);
  eff_mm->SetMarkerSize(1.5);
  eff_mm->SetMarkerColor(kBlack);

  TLegend* legend = new TLegend( 0.2, 0.2, 0.6, 0.4 );
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  legend->SetFillColor(0);
  legend->AddEntry( eff_ee, "Electrons (cut = 0.1)", "P" );
  legend->AddEntry( eff_mm, "Muons (cut = 0.2)", "P" );
  legend->Draw("same");

  eff_ee->Draw("p same");
  eff_mm->Draw("p same");

  ZGDrawTools::addLabels( c1, -1., "CMS Simulation" );

  gPad->RedrawAxis();

  c1->SaveAs("eff_miniIso.eps");
  c1->SaveAs("eff_miniIso.pdf");

  return 0;

}
