#include "TFile.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "../interface/ZGDrawTools.h"



int main() {

  ZGDrawTools::setStyle();


  TFile* file = TFile::Open("genSignalShapes/genSignalShapes.root");


  TGraphErrors* gr_mean_1p4 = (TGraphErrors*)file->Get("bw_mean_1p4");
  TGraphErrors* gr_mean_5p6 = (TGraphErrors*)file->Get("bw_mean_5p6");

  TGraphErrors* gr_width_1p4 = (TGraphErrors*)file->Get("bw_width_1p4");
  TGraphErrors* gr_width_5p6 = (TGraphErrors*)file->Get("bw_width_5p6");

  TCanvas* c1 = new TCanvas( "c1", "c1", 600, 600 );
  c1->cd();

  TH2D* h2_axes = new TH2D("axes", "", 10, 300., 1000., 10, 300., 1000.);
  h2_axes->SetXTitle( "Generated Mass [GeV]" );
  h2_axes->SetYTitle( "Fitted Mass [GeV]" );
  h2_axes->Draw();

  TLine* diag = new TLine( 300., 300., 1000., 1000.);
  diag->SetLineColor(46);
  diag->SetLineWidth(2);
  diag->Draw("same");

  gr_mean_1p4->SetMarkerStyle(24);
  gr_mean_1p4->SetMarkerSize(1.3);

  gr_mean_5p6->SetMarkerStyle(20);
  gr_mean_5p6->SetMarkerSize(1.3);

  gr_mean_1p4->Draw("psame");
  gr_mean_5p6->Draw("psame");

  TLegend* legend = new TLegend(0.2, 0.7, 0.5, 0.9);
  legend->SetFillColor(0);
  legend->SetTextSize(0.035);
  legend->SetTextFont(42);
  legend->AddEntry( diag, "Expected", "L" );
  legend->AddEntry( gr_mean_1p4, "W = 1.4%", "P" );
  legend->AddEntry( gr_mean_5p6, "W = 5.6%", "P" );
  legend->Draw("same");

  gPad->RedrawAxis();

  ZGDrawTools::addLabels( c1, -1, "CMS Simulation" );

  c1->SaveAs( "genSignalShapes/mean.eps" );
  c1->SaveAs( "genSignalShapes/mean.pdf" );

  c1->Clear();


  TH2D* h2_axes2 = new TH2D("axes2", "", 10, 300., 1000., 10, 0., 0.1);
  h2_axes2->SetXTitle( "Generated #Gamma/M" );
  h2_axes2->SetYTitle( "Fitted #Gamma/M" );
  h2_axes2->Draw();

  TLine* line_1p4 = new TLine( 300., 0.014, 1000., 0.014);
  line_1p4->SetLineColor(46);
  line_1p4->SetLineWidth(2);
  line_1p4->Draw("same");

  TLine* line_5p6 = new TLine( 300., 0.056, 1000., 0.056);
  line_5p6->SetLineColor(46);
  line_5p6->SetLineWidth(2);
  line_5p6->Draw("same");

  gr_width_1p4->SetMarkerStyle(24);
  gr_width_1p4->SetMarkerSize(1.3);

  gr_width_5p6->SetMarkerStyle(20);
  gr_width_5p6->SetMarkerSize(1.3);

  gr_width_1p4->Draw("psame");
  gr_width_5p6->Draw("psame");

  legend->Draw("same");

  gPad->RedrawAxis();

  ZGDrawTools::addLabels( c1, -1, "CMS Simulation" );

  c1->SaveAs( "genSignalShapes/width.eps" );
  c1->SaveAs( "genSignalShapes/width.pdf" );

  

  return 0;

}
