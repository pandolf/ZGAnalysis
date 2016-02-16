#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLegend.h"


#include "../interface/ZGDrawTools.h"



int main() {

  ZGDrawTools::setStyle();

  TFile* file_74 = TFile::Open("EventYields_v0_eth74X/trees.root");
  TFile* file_76 = TFile::Open("EventYields_v0_eth76X/trees.root");

  TTree* tree_74 = (TTree*)file_74->Get("data");
  TTree* tree_76 = (TTree*)file_76->Get("data");

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();
 

  int nBins = 15;
  float xMin = 200.;
  float xMax = 500.;
  float binWidth = (xMax-xMin)/((float)nBins);

  TH1D* h1_74 = new TH1D( "h1_74", "", nBins, xMin, xMax );
  TH1D* h1_76 = new TH1D( "h1_76", "", nBins, xMin, xMax );

  h1_74->SetXTitle( "M(Z#gamma) [GeV]");
  h1_74->SetYTitle( Form("Events / (%.0f GeV)", binWidth) );

  h1_74->SetLineWidth(2);
  h1_76->SetLineWidth(2);

  tree_74->Project("h1_74", "boss_mass");
  tree_76->Project("h1_76", "boss_mass");
  h1_74->Draw("same");
  h1_76->SetLineColor(kRed);
  h1_76->Draw("same");

  ZGDrawTools::addLabels( c1, 2.6 );

  TLegend* legend = new TLegend( 0.55, 0.75, 0.9, 0.9 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  legend->AddEntry( h1_74, "74X Data", "L" );
  legend->AddEntry( h1_76, "76X Data", "L" );
  legend->Draw("same");

  gPad->RedrawAxis();

  c1->SaveAs( "mass_74vs76.eps" );
  c1->SaveAs( "mass_74vs76.pdf" );

  return 0;

}
