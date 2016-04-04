#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TLegend.h"

#include "../interface/ZGDrawTools.h"




void computeHLTEfficiency( TFile* file, int mass, TGraph* gr_eff_ee, TGraph* gr_eff_mm );


int main() {

  ZGDrawTools::setStyle();

  TFile* file = TFile::Open("EventYields_v0_eth74X_noHLT/trees.root");
  std::vector<int> masses;
  masses.push_back(350);
  masses.push_back(450);
  masses.push_back(600);
  masses.push_back(750);
  masses.push_back(900);
  masses.push_back(1500);
  masses.push_back(2000);


  TGraph* gr_eff_mm = new TGraph(0);
  TGraph* gr_eff_ee = new TGraph(0);

  for( unsigned i=0; i<masses.size(); ++i ) {

    computeHLTEfficiency( file, masses[i], gr_eff_ee, gr_eff_mm );

  }

  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();

  TH2D* h2_axes = new TH2D("axes", "", 10, 200., 2100., 10, 0.5, 1.1 );
  h2_axes->SetXTitle( "Generated Signal Mass [GeV]");
  h2_axes->SetYTitle( "Trigger Efficiency" );
  h2_axes->Draw();

  TLine* line_one = new TLine( 200., 1., 2100., 1.);
  line_one->Draw("same");

  gr_eff_mm->SetMarkerStyle(20);
  gr_eff_mm->SetMarkerSize(1.5);
  gr_eff_mm->SetMarkerColor(38);

  gr_eff_ee->SetMarkerStyle(20);
  gr_eff_ee->SetMarkerSize(1.5);
  gr_eff_ee->SetMarkerColor(46);

  gr_eff_ee->Draw("p same");
  gr_eff_mm->Draw("p same");

  TLegend* legend = new TLegend( 0.2, 0.2, 0.5, 0.4 );
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  legend->SetFillColor(0);
  legend->AddEntry( gr_eff_ee, "Electron Channel", "P");
  legend->AddEntry( gr_eff_mm, "Muon Channel", "P");
  legend->Draw("same");

  ZGDrawTools::addLabels(c1, -1., "CMS Simulation");

  gPad->RedrawAxis();

  c1->SaveAs( "EventYields_v0_eth74X_noHLT/signalHLTEfficiency.eps" );
  c1->SaveAs( "EventYields_v0_eth74X_noHLT/signalHLTEfficiency.pdf" );

  return 0;

}



void computeHLTEfficiency( TFile* file, int mass, TGraph* gr_eff_ee, TGraph* gr_eff_mm ) {

  TTree* tree = (TTree*)file->Get(Form("XZg_Spin0_ZToLL_W_0p014_M_%d", mass) );

  int denom_mm = tree->GetEntries("leptType==13");
  int num_mm   = tree->GetEntries("leptType==13 && passHLT");

  int denom_ee = tree->GetEntries("leptType==11");
  int num_ee   = tree->GetEntries("leptType==11 && passHLT");
  
  float eff_ee = (float)num_ee/(float)denom_ee;
  float eff_mm = (float)num_mm/(float)denom_mm;

  gr_eff_ee->SetPoint( gr_eff_ee->GetN(), mass, eff_ee );
  gr_eff_mm->SetPoint( gr_eff_mm->GetN(), mass, eff_mm );

}
