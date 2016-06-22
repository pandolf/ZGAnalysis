#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLegend.h"


#include "../interface/ZGDrawTools.h"



int main( int argc, char* argv[] ) {

  ZGDrawTools::setStyle();

  if( argc==1 ) {
    std::cout << "USAGE: ./drawOptimizeCut [lumi]" << std::endl;
    exit(1);
  }

  std::string dir(argv[1]);
  int lumi = atoi(lumi_str.c_str());
  std::string dir(Form("optCuts_%sfb", lumi_str.c_str()));

  TFile* file = TFile::Open( Form("%s/file_opt.root", dir.c_str()) );
  
  
  TH1D* sig_350 = (TH1D*)file->Get("sig_350");
  TH1D* sig_450 = (TH1D*)file->Get("sig_450");
  TH1D* sig_600 = (TH1D*)file->Get("sig_600");
  TH1D* sig_750 = (TH1D*)file->Get("sig_750");

  sig_350->SetMarkerColor(46);
  sig_450->SetMarkerColor(29);
  sig_600->SetMarkerColor(kBlack);
  sig_750->SetMarkerColor(38);

  sig_350->SetMarkerSize(1.3);
  sig_450->SetMarkerSize(1.3);
  sig_600->SetMarkerSize(1.3);
  sig_750->SetMarkerSize(1.3);

  sig_350->SetMarkerStyle(20);
  sig_450->SetMarkerStyle(21);
  sig_600->SetMarkerStyle(24);
  sig_750->SetMarkerStyle(25);


  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();

  float yMax350 = sig_350->GetMaximum();
  float yMax750 = sig_750->GetMaximum();
  float yMax = (yMax350>yMax750) ? yMax350 : yMax750;
  
  TH2D* h2_axes = new TH2D( "axes", "", 10, 0., 1., 10, 0., yMax*1.4 );
  h2_axes->SetXTitle( "Threshold on p_{T}(#gamma)/M(Z#gamma)" );
  h2_axes->SetYTitle( "s / #sqrt{b} " );
  h2_axes->Draw("same");

  TLegend* legend = new TLegend( 0.6, 0.6, 0.9, 0.9 );
  legend->SetFillColor(0); 
  legend->SetTextSize(0.035);
  legend->SetHeader( Form("L = %d fb^{-1}", lumi) );
  legend->AddEntry( sig_350, "M = 350 GeV", "P" );
  legend->AddEntry( sig_450, "M = 450 GeV", "P" );
  legend->AddEntry( sig_600, "M = 600 GeV", "P" );
  legend->AddEntry( sig_750, "M = 750 GeV", "P" );
  legend->Draw("same");

  sig_350->Draw("Psame");
  sig_450->Draw("Psame");
  sig_600->Draw("Psame");
  sig_750->Draw("Psame");

  ZGDrawTools::addLabels( c1, -1., "CMS Simulation");

  gPad->RedrawAxis();

  c1->SaveAs( Form( "%s/optPlot.eps", dir.c_str()) );
  c1->SaveAs( Form( "%s/optPlot.pdf", dir.c_str()) );


  return 0;

}


