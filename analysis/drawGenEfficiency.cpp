#include <cmath>

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


float xMin = 300.;
float xMax = 1000.;



void drawEfficiency( const std::string& outdir, TFile* file );
void drawVsMass( const std::string& outdir, TFile* file );
void drawSingleGraph( const std::string& outdir, TGraphErrors* graph, float yMin, float yMax, const std::string& axisName, float lineY=-999. );



int main() {

  ZGDrawTools::setStyle();

  TFile* file = TFile::Open("provaGenEfficiency.root");

  std::string outdir = "plotsGenEfficiency";
  system( Form("mkdir -p %s", outdir.c_str()) );

  drawEfficiency( outdir, file );  

  drawVsMass( outdir, file );

  return 0;

}


void drawEfficiency( const std::string& outdir, TFile* file ) {

  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();

  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., 1.001 );
  h2_axes->SetXTitle("Generated Z#gamma Mass [GeV]");
  h2_axes->SetYTitle("Reconstruction Efficiency");
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
  h1_band->SetFillColor(18);
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



void drawVsMass( const std::string& outdir, TFile* file ) {

  TH1D* h1_all = (TH1D*)file->Get("massReso");

  TGraphErrors* gr_resp = new TGraphErrors(0);
  TGraphErrors* gr_reso = new TGraphErrors(0);
  gr_resp->SetName( "resp" );
  gr_reso->SetName( "reso" );

  for( unsigned i=0; i<100; i++ ) {

    TCanvas* c1 = new TCanvas( Form("fit_%d", i), "", 600, 600 );
    c1->cd();

    TH1D* h1_reso = (TH1D*)file->Get( Form("reso_%d", i) );
    if( h1_reso==0 ) break;

    TF1* f1 = new TF1("gaussian", "gaus", -0.5, 0.5);
    f1->SetParameters( h1_reso->Integral(), h1_reso->GetMean(), h1_reso->GetRMS() );
    h1_reso->Fit( f1, "LQR" );

    h1_reso->Draw();
    c1->SaveAs( Form("%s/fit_%d.eps", outdir.c_str(), i) );
    c1->SaveAs( Form("%s/fit_%d.pdf", outdir.c_str(), i) );

    float xmass = h1_all->GetBinCenter( i+1 );
    float mean = f1->GetParameter(1);
    float mean_err = f1->GetParError(1);
    gr_resp->SetPoint(i, xmass, mean );
    gr_resp->SetPointError(i, h1_all->GetBinWidth( i+1 )/sqrt(12.), mean_err );

    float sigma = f1->GetParameter(2);
    float sigma_err = f1->GetParError(2);
    gr_reso->SetPoint(i, xmass, sigma );
    gr_reso->SetPointError(i, h1_all->GetBinWidth( i+1 )/sqrt(12.), sigma_err );
    //gr_reso->SetPointError(i, h1_all->GetBinWidth( i+1 )/sqrt(12.), sqrt( sigma_err*sigma_err/(mean*mean) + sigma*sigma*mean_err*mean_err/(mean*mean*mean*mean) ) );

    delete c1;

  }

  drawSingleGraph( outdir, gr_resp, -0.1, 0.1, "(M_{gen} - M_{reco}) / M_{gen}", 0. );
  drawSingleGraph( outdir, gr_reso, 0., 0.2, "Mass Resolution" );

}


void drawSingleGraph( const std::string& outdir, TGraphErrors* graph, float yMin, float yMax, const std::string& axisName, float lineY ) {


  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();

  TH2D* h2_axes = new TH2D( "axes", "", 10, xMin, xMax, 10, yMin, yMax );
  h2_axes->SetXTitle( "Generated Z#gamma Mass [GeV]" );
  h2_axes->SetYTitle( axisName.c_str() );
  h2_axes->Draw();

  TLine* line = new TLine( xMin, lineY, xMax, lineY );
  line->Draw("same");
  

  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(2.);
  graph->Draw("p same" );
  
  ZGDrawTools::addLabels( c1, -1., "CMS Simulation" );

  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/%s.eps", outdir.c_str(), graph->GetName()) );
  c1->SaveAs( Form("%s/%s.pdf", outdir.c_str(), graph->GetName()) );

  delete c1;
  delete h2_axes;

}
