#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLegend.h"


#include "../interface/ZGDrawTools.h"


void drawVar( TTree* tree_74, TTree* tree_76, const std::string& savename, const std::string& name, int nBins, float xMin, float xMax, const std::string& selection="" );

int main() {

  ZGDrawTools::setStyle();

  TFile* file_74 = TFile::Open("EventYields_v0_eth74X/trees.root");
  TFile* file_76 = TFile::Open("EventYields_v0_eth76X_2/trees.root");

  TTree* tree_74 = (TTree*)file_74->Get("data");
  TTree* tree_76 = (TTree*)file_76->Get("data");

  drawVar( tree_74, tree_76, "boss_mass", "boss_mass", 15, 200., 500.);
  drawVar( tree_74, tree_76, "gamma_pt","gamma_pt", 25, 40., 300.);
  drawVar( tree_74, tree_76, "mZee", "z_mass", 80, 50., 130., "leptType==11");
  drawVar( tree_74, tree_76, "mZmm", "z_mass", 80, 50., 130., "leptType==13");

  return 0;

}

void drawVar( TTree* tree_74, TTree* tree_76, const std::string& savename, const std::string& name, int nBins, float xMin, float xMax, const std::string& selection ) {

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();
 

  float binWidth = (xMax-xMin)/((float)nBins);

  TH1D* h1_74 = new TH1D( "h1_74", "", nBins, xMin, xMax );
  TH1D* h1_76 = new TH1D( "h1_76", "", nBins, xMin, xMax );

  h1_76->SetXTitle( name.c_str() );
  //h1_76->SetXTitle( "M(Z#gamma) [GeV]");
  h1_76->SetYTitle( "Events" );
  //h1_76->SetYTitle( Form("Events / (%.0f GeV)", binWidth) );

  h1_74->SetLineWidth(2);
  h1_76->SetLineWidth(2);

  tree_74->Project("h1_74", name.c_str(), selection.c_str());
  tree_76->Project("h1_76", name.c_str(), selection.c_str());
  h1_76->SetLineColor(kRed);
  h1_76->Draw("same");
  h1_74->Draw("same");

  ZGDrawTools::addLabels( c1, 2.6 );

  TLegend* legend = new TLegend( 0.65, 0.75, 0.9, 0.9 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  legend->AddEntry( h1_74, "74X Data", "L" );
  legend->AddEntry( h1_76, "76X Data", "L" );
  legend->Draw("same");

  gPad->RedrawAxis();

  c1->SaveAs( Form("%s_74vs76.eps", savename.c_str()) );
  c1->SaveAs( Form("%s_74vs76.pdf", savename.c_str()) );

  c1->SetLogy();

  c1->SaveAs( Form("%s_74vs76_log.eps", savename.c_str()) );
  c1->SaveAs( Form("%s_74vs76_log.pdf", savename.c_str()) );
  
  delete c1;
  delete h1_74;
  delete h1_76;

}
