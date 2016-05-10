#include <fstream>

#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphAsymmErrors.h"

#include "../interface/ZGDrawTools.h"
#include "../interface/ZGConfig.h"


void scaleGraph( TGraphAsymmErrors* gr );


int main( int argc, char* argv[] ) {


  ZGDrawTools::setStyle();


  std::string axisName = "95\% CL UL on #sigma #times BR(A#rightarrowZ#gamma) [fb]";

  TFile* file = TFile::Open("limit_plots.root");
  TGraphAsymmErrors* gr_exp        = (TGraphAsymmErrors*)file->Get("expected");
  TGraphAsymmErrors* gr_exp_1sigma = (TGraphAsymmErrors*)file->Get("oneSigma");
  TGraphAsymmErrors* gr_exp_2sigma = (TGraphAsymmErrors*)file->Get("twoSigma");
  TGraphAsymmErrors* gr_obs        = (TGraphAsymmErrors*)file->Get("observed");
  
  scaleGraph( gr_exp );
  scaleGraph( gr_exp_1sigma );
  scaleGraph( gr_exp_2sigma );
  scaleGraph( gr_obs        );

  gr_obs->SetLineWidth(2);

  gr_exp->SetLineWidth(2);
  gr_exp->SetLineStyle(2);

  gr_exp_1sigma->SetLineWidth(0);
  gr_exp_1sigma->SetFillColor(8);

  gr_exp_2sigma->SetLineWidth(0);
  gr_exp_2sigma->SetFillColor(219);

  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();

  float yMax = 150.;

  TH2D* h2_axes = new TH2D("axes", "", 10, 200., 1200., 10, 0., yMax );
  h2_axes->SetYTitle( axisName.c_str() );
  h2_axes->SetXTitle( "Resonance Mass [GeV]");
  h2_axes->GetXaxis()->SetNdivisions(1006);
  h2_axes->Draw();


  gr_exp_2sigma->Draw("E3 same");
  gr_exp_1sigma->Draw("E3 same");
  gr_exp       ->Draw("L  same");
  gr_obs       ->Draw("L  same");


  gr_exp_1sigma->SetLineWidth(2);
  gr_exp_1sigma->SetLineStyle(2);
  gr_exp_2sigma->SetLineWidth(2);
  gr_exp_2sigma->SetLineStyle(2);

  TLegend* legend;
  std::string title = "Narrow Signal Model";
  legend = new TLegend( 0.55, 0.6 , 0.9, 0.9 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  legend->SetHeader(title.c_str());
  legend->AddEntry( gr_obs, "Observed", "L" );
  legend->AddEntry( gr_exp_1sigma, "Expected #pm 1#sigma", "LF" );
  legend->AddEntry( gr_exp_2sigma, "Expected #pm 2#sigma", "LF" );
  legend->Draw("same");

  TPaveText* label_cms = new TPaveText(0.143,0.96,0.27,0.965, "brNDC");
  label_cms->SetBorderSize(0);
  label_cms->SetFillColor(kWhite);
  label_cms->SetTextSize(0.042);
  label_cms->SetTextAlign(11); // align left
  label_cms->SetTextFont(42);  
  label_cms->AddText( "CMS" );
  label_cms->Draw("same");

  TPaveText* label_top = new TPaveText(0.4,0.959,0.975,0.963, "brNDC");
  label_top->SetBorderSize(0);
  label_top->SetFillColor(kWhite);
  label_top->SetTextSize(0.038);
  label_top->SetTextAlign(31); // align right
  label_top->SetTextFont(42);  // label_top->SetTextFont(62);
  label_top->AddText("19.7 fb^{-1} (8 TeV)");
  label_top->Draw("same");


  gPad->RedrawAxis();

  c1->SaveAs( "limitPlot8TeV.eps" );
  c1->SaveAs( "limitPlot8TeV.pdf" );

  return 0;

}



void scaleGraph( TGraphAsymmErrors* gr ) {

  for( int i=0; i<gr->GetN(); ++i ) {
    double x,y;
    gr->GetPoint(i,x,y);
    double yerrl = gr->GetErrorYlow(i);
    double yerrh = gr->GetErrorYhigh(i);
    gr->SetPoint(i,x,y/0.03366);
    gr->SetPointEYlow (i, yerrl/0.03366);
    gr->SetPointEYhigh(i, yerrh/0.03366);
  }

}
