#include "TFile.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "../interface/ZGDrawTools.h"


void compareParameter( TFile* file, TFile* outfile, const std::string& name, const std::string& axisName );


int main() {

  ZGDrawTools::setStyle();


  TFile* file = TFile::Open("genSignalShapes/genSignalShapes.root");
  TFile* outfile = TFile::Open("genSignalShapes/genSignalShapesParameters.root", "recreate");

  compareParameter( file, outfile, "mean" , "Gaussian #mu [GeV]" );
  compareParameter( file, outfile, "sigma", "Gaussian #sigma [GeV]" );
  compareParameter( file, outfile, "width", "Gaussian #sigma/#mu" );
  compareParameter( file, outfile, "alpha1", "CB #alpha_{1}" );
  compareParameter( file, outfile, "alpha2", "CB #alpha_{2}" );
  compareParameter( file, outfile, "n1", "CB n_{1}" );
  compareParameter( file, outfile, "n2", "CB n_{2}" );

  outfile->Close();

  return 0;

}


void compareParameter( TFile* file, TFile* outfile, const std::string& name, const std::string& axisName ) {

  TGraphErrors* gr_1p4 = (TGraphErrors*)file->Get(Form("%s_1p4", name.c_str()) );
  TGraphErrors* gr_5p6 = (TGraphErrors*)file->Get(Form("%s_5p6", name.c_str()) );


  TCanvas* c1 = new TCanvas( "c1", "c1", 600, 600 );
  c1->cd();

  float yMax = 0.;
  for( int i=0; i<gr_5p6->GetN(); ++i ) {
    Double_t x, y;
    gr_5p6->GetPoint(i, x, y);
    if( y>yMax ) yMax = y;
  }

  float xMin = 350.;
  float xMax = 1000.;

  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., 1.3*yMax );
  h2_axes->SetXTitle( "Generated Mass [GeV]" );
  h2_axes->SetYTitle( axisName.c_str() );
  h2_axes->Draw();

  TF1* f1_1p4 = new TF1( Form("f1_%s", gr_1p4->GetName()), "[0] + [1]*x", xMin, xMax );
  f1_1p4->SetLineColor(38);

  gr_1p4->Fit( f1_1p4, "QR" );
  TH1D* band_1p4 = ZGDrawTools::getBand(f1_1p4);

  gr_1p4->SetMarkerStyle(20);
  gr_1p4->SetMarkerSize(2);
  gr_1p4->SetMarkerColor(38);
  gr_1p4->SetLineColor(38);
  gr_1p4->SetLineWidth(2);


  TF1* f1_5p6 = new TF1( Form("f1_%s", gr_5p6->GetName()), "[0] + [1]*x", xMin, xMax );
  f1_5p6->SetLineColor(46);

  gr_5p6->Fit( f1_5p6, "QR" );
  TH1D* band_5p6 = ZGDrawTools::getBand(f1_5p6);

  gr_5p6->SetMarkerStyle(20);
  gr_5p6->SetMarkerSize(2);
  gr_5p6->SetMarkerColor(46);
  gr_5p6->SetLineColor(46);
  gr_5p6->SetLineWidth(2);

  TLegend* legend = new TLegend(0.2, 0.75, 0.5, 0.9);
  legend->SetFillColor(0);
  legend->SetTextSize(0.035);
  legend->SetTextFont(42);
  legend->AddEntry( gr_1p4, "W = 1.4%", "PL" );
  legend->AddEntry( gr_5p6, "W = 5.6%", "PL" );
  legend->Draw("same");


  band_1p4->Draw("c e3 same");
  band_5p6->Draw("c e3 same");

  gr_1p4->Draw("psame");
  gr_5p6->Draw("psame");

  gPad->RedrawAxis();

  ZGDrawTools::addLabels( c1, -1, "CMS Simulation" );

  c1->SaveAs( Form("genSignalShapes/%s.eps", name.c_str()) );
  c1->SaveAs( Form("genSignalShapes/%s.pdf", name.c_str()) );

  outfile->cd();
  f1_1p4->Write();
  f1_5p6->Write();

  delete c1;
  delete h2_axes;

}
