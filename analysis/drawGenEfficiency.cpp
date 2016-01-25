#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TEfficiency.h"
#include "TVirtualFitter.h"


#include "../interface/ZGDrawTools.h"


void drawEfficiency( const std::string& outdir, TFile* file );



int main() {

  ZGDrawTools::setStyle();

  TFile* file = TFile::Open("provaGenEfficiency.root");

  std::string outdir = "plotsGenEfficiency";
  system( Form("mkdir -p %s", outdir.c_str()) );

  drawEfficiency( outdir, file );  

  return 0;

}


void drawEfficiency( const std::string& outdir, TFile* file ) {

  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();

  float xMin = 300.;
  float xMax = 1000.;

  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., 1.001 );
  h2_axes->SetXTitle("Generated Z#gamma Mass [GeV]");
  h2_axes->SetYTitle("Efficiency");
  h2_axes->Draw();

  TEfficiency* eff = (TEfficiency*)file->Get("eff");
  TGraphAsymmErrors* gr_eff = eff->CreateGraph();
  gr_eff->SetMarkerStyle(20);
  gr_eff->SetMarkerSize(1.3);

  TF1* line = new TF1( "line", "[0] + [1]*x", xMin, xMax );
  line->SetLineColor(46);
  line->SetLineWidth(2);
  gr_eff->Fit( line, "QR0" );

  TH1D* h1_band = new TH1D("band", "", 500, 300., 1000.);
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(h1_band, 0.68);
  h1_band->SetFillColor(42);
  h1_band->Draw("C E3 same");
  line->Draw("same");

  gr_eff->Draw("p same");

  

  ZGDrawTools::addLabels( c1, -1., "CMS Simulation");

  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/eff.eps", outdir.c_str()) );
  c1->SaveAs( Form("%s/eff.pdf", outdir.c_str()) );

  delete c1;
  delete h2_axes;

}
