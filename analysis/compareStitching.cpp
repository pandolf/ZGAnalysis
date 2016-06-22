#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLegend.h"


#include "../interface/ZGDrawTools.h"


void drawVar( const std::string& outdir, const std::string& saveName, const std::string& varName, const std::string& axisName, int nBins, float xMin, float xMax, TTree* tree_old, TTree* tree_stitch, TTree* tree_stitch_ptWeight, const std::string& selection="" );

int main() {

  ZGDrawTools::setStyle();

  TFile* file_old             = TFile::Open("$DCAP/pnfs/psi.ch/cms/trivcat/store/user/pandolf/ZGproduction/80X/PostProcessed/mc2016_v1_pu/skimAndPrune/ZGTo2LG_post_skim.root");
  TFile* file_stitch          = TFile::Open("ZGTo2LG_post_skim_stitch.root");
  TFile* file_stitch_ptWeight = TFile::Open("ZGTo2LG_post_skim_stitch_ptWeight.root");

  TTree* tree_old             = (TTree*)file_old            ->Get("mt2");
  TTree* tree_stitch          = (TTree*)file_stitch         ->Get("mt2");
  TTree* tree_stitch_ptWeight = (TTree*)file_stitch_ptWeight->Get("mt2");


  std::string outdir = "compareStitch";
  system( Form("mkdir -p %s", outdir.c_str()));

  drawVar( outdir, "gamma_pt", "gamma_pt[0]", "Photon p_{T} [GeV]", 25, 80., 380., tree_old, tree_stitch, tree_stitch_ptWeight, "nlep==2 && ngamma>0" );
  drawVar( outdir, "lept0_pt", "lep_pt[0]", "Leading Lepton p_{T} [GeV]", 25, 0., 300., tree_old, tree_stitch, tree_stitch_ptWeight, "nlep==2 && ngamma>0" );

  return 0;

}


void drawVar( const std::string& outdir, const std::string& saveName, const std::string& varName, const std::string& axisName, int nBins, float xMin, float xMax, TTree* tree_old, TTree* tree_stitch, TTree* tree_stitch_ptWeight, const std::string& selection ) {


  TH1D* h1_old             = new TH1D(Form("%s_old"            , varName.c_str()), "", nBins, xMin, xMax );
  TH1D* h1_stitch          = new TH1D(Form("%s_stitch"         , varName.c_str()), "", nBins, xMin, xMax );
  TH1D* h1_stitch_ptWeight = new TH1D(Form("%s_stitch_ptWeight", varName.c_str()), "", nBins, xMin, xMax );

  h1_old            ->Sumw2();
  h1_stitch         ->Sumw2();
  h1_stitch_ptWeight->Sumw2();

  tree_old            ->Project( h1_old            ->GetName(), varName.c_str(), Form("evt_scale1fb*(%s)", selection.c_str()) );
  tree_stitch         ->Project( h1_stitch         ->GetName(), varName.c_str(), Form("evt_scale1fb*(%s)", selection.c_str()) );
  tree_stitch_ptWeight->Project( h1_stitch_ptWeight->GetName(), varName.c_str(), Form("evt_scale1fb*(%s)", selection.c_str()) );


  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->SetLogy();
  c1->cd();

  float yMax = h1_old->GetMaximum()*1.3;

  TH2D* h2_axes = new TH2D( "axes", "", 10, xMin, xMax, 10, 0.07, yMax );
  h2_axes->SetXTitle( axisName.c_str() );
  h2_axes->SetYTitle( "Events" );
  h2_axes->Draw();

  h1_old            ->SetLineColor(kBlack);
  h1_stitch         ->SetLineColor(kBlue);
  h1_stitch_ptWeight->SetLineColor(kRed);

  h1_old            ->SetLineWidth(2.);
  h1_stitch         ->SetLineWidth(2.);
  h1_stitch_ptWeight->SetLineWidth(2.);

  TLegend* legend = new TLegend( 0.55, 0.7, 0.9, 0.9 );
  legend->SetFillColor(0);
  legend->SetTextSize( 0.035 );
  legend->AddEntry( h1_old, "Inclusive Sample", "L" );
  legend->AddEntry( h1_stitch, "Stitched", "L" );
  legend->AddEntry( h1_stitch_ptWeight, "+ p_{T} k-factor", "L" );
  legend->Draw("same");

  h1_old            ->Draw("same");
  h1_stitch         ->Draw("same");
  h1_stitch_ptWeight->Draw("same");

  gPad->RedrawAxis();

  ZGDrawTools::addLabels( c1, -1., "CMS Simulation" );

  c1->SaveAs( Form("%s/%s.eps", outdir.c_str(), saveName.c_str()) );
  c1->SaveAs( Form("%s/%s.pdf", outdir.c_str(), saveName.c_str()) );

  delete c1;
  delete h2_axes;
  delete h1_old;
  delete h1_stitch;
  delete h1_stitch_ptWeight;

}

