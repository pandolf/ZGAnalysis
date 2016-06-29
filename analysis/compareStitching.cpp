#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLegend.h"


#include "../interface/ZGDrawTools.h"


void drawVar( const std::string& outdir, const std::string& saveName, const std::string& varName, const std::string& axisName, int nBins, float xMin, float xMax, TTree* tree_old, TTree* tree_stitch, const std::string& selection="" );

int main() {

  ZGDrawTools::setStyle();

  TFile* file_old    = TFile::Open("EventYields_v0_2016v3/trees.root");
  TFile* file_stitch = TFile::Open("EventYields_v0_stitch/trees.root");

  TTree* tree_old             = (TTree*)file_old            ->Get("zg");
  TTree* tree_stitch          = (TTree*)file_stitch         ->Get("zg");


  std::string outdir = "compareStitch";
  system( Form("mkdir -p %s", outdir.c_str()));

  drawVar( outdir, "gamma_pt", "gamma_pt" , "Photon p_{T} [GeV]", 50, 80., 880.,        tree_old, tree_stitch );
  drawVar( outdir, "lept0_pt", "lept0_pt" , "Leading Lepton p_{T} [GeV]", 50, 0., 600., tree_old, tree_stitch );
  drawVar( outdir, "mZg"     , "boss_mass", "M(Z#gamma) [GeV]", 50, 200., 2200.,        tree_old, tree_stitch );
  //drawVar( outdir, "mZg"     , "boss_mass", "M(Z#gamma) [GeV]", 50, 200., 2200., tree_old, tree_stitch, tree_stitch_ptWeight );

  //drawVar( outdir, "gamma_pt", "gamma_pt[0]", "Photon p_{T} [GeV]", 25, 80., 380., tree_old, tree_stitch, tree_stitch_ptWeight, "nlep==2 && ngamma>0" );
  //drawVar( outdir, "lept0_pt", "lep_pt[0]", "Leading Lepton p_{T} [GeV]", 25, 0., 300., tree_old, tree_stitch, tree_stitch_ptWeight, "nlep==2 && ngamma>0" );

  return 0;

}


void drawVar( const std::string& outdir, const std::string& saveName, const std::string& varName, const std::string& axisName, int nBins, float xMin, float xMax, TTree* tree_old, TTree* tree_stitch, const std::string& selection ) {


  TH1D* h1_old             = new TH1D(Form("%s_old"            , varName.c_str()), "", nBins, xMin, xMax );
  TH1D* h1_stitch          = new TH1D(Form("%s_stitch"         , varName.c_str()), "", nBins, xMin, xMax );

  h1_old            ->Sumw2();
  h1_stitch         ->Sumw2();

  if( selection=="" ) {
    tree_old            ->Project( h1_old            ->GetName(), varName.c_str(), Form("weight") );
    tree_stitch         ->Project( h1_stitch         ->GetName(), varName.c_str(), Form("weight") );
  } else {
    tree_old            ->Project( h1_old            ->GetName(), varName.c_str(), Form("weight*(%s)", selection.c_str()) );
    tree_stitch         ->Project( h1_stitch         ->GetName(), varName.c_str(), Form("weight*(%s)", selection.c_str()) );
  }


  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->SetLogy();
  c1->cd();

  float yMax = h1_old->GetMaximum()*1.3;

  TH2D* h2_axes = new TH2D( "axes", "", 10, xMin, xMax, 10, 0.0007, yMax );
  h2_axes->SetXTitle( axisName.c_str() );
  h2_axes->SetYTitle( "Events" );
  h2_axes->Draw();

  h1_old            ->SetLineColor(kBlack);
  h1_stitch         ->SetLineColor(kRed);

  h1_old            ->SetLineWidth(2.);
  h1_stitch         ->SetLineWidth(2.);

  TLegend* legend = new TLegend( 0.55, 0.75, 0.9, 0.9 );
  legend->SetFillColor(0);
  legend->SetTextSize( 0.035 );
  legend->AddEntry( h1_old, "Inclusive Sample", "L" );
  legend->AddEntry( h1_stitch, "Stitched", "L" );
  legend->Draw("same");

  h1_old            ->Draw("same");
  h1_stitch         ->Draw("same");

  gPad->RedrawAxis();

  ZGDrawTools::addLabels( c1, -1., "CMS Simulation" );

  c1->SaveAs( Form("%s/%s_stitch.eps", outdir.c_str(), saveName.c_str()) );
  c1->SaveAs( Form("%s/%s_stitch.pdf", outdir.c_str(), saveName.c_str()) );

  delete c1;
  delete h2_axes;
  delete h1_old;
  delete h1_stitch;

}

