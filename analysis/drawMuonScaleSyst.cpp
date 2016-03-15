#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "../interface/ZGDrawTools.h"
#include "../interface/ZGConfig.h"



int main( int argc, char* argv[] ) {

  if( argc<2 ) {
    std::cout << "USAGE: ./drawMuonScaleSyst [configFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  ZGDrawTools::setStyle();

  std::string configFileName(argv[1]);
  ZGConfig cfg(configFileName);

  std::string dir = cfg.getEventYieldDir() + "/muonScaleSyst";
  

  TFile* file = TFile::Open( Form("%s/muonScaleSyst.root", dir.c_str()) );

  TGraph* gr_muScale  = (TGraph*)file->Get("gr_muScaleSyst" );
  TGraph* gr_egmScale = (TGraph*)file->Get("gr_egmScaleSyst");

  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();

  TH2D* h2_axes = new TH2D("axes", "", 10, 200., 2100., 10, 0., 0.07);
  h2_axes->SetXTitle( "Generated Signal Mass [GeV]" );
  h2_axes->SetYTitle( "Muon Scale Systematic Uncertainty" );
  h2_axes->Draw();

  gr_muScale->SetMarkerStyle(20); 
  gr_muScale->SetMarkerSize(1.5);
  gr_muScale->SetMarkerColor(kBlack);

  gr_muScale->Draw("psame");

  ZGDrawTools::addLabels(c1, -1., "CMS Simulation");

  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/muonScaleSyst.eps", dir.c_str()) );
  c1->SaveAs( Form("%s/muonScaleSyst.pdf", dir.c_str()) );

  c1->Clear();


  h2_axes->Draw();

  gr_egmScale->SetMarkerStyle(20); 
  gr_egmScale->SetMarkerSize(1.5);
  gr_egmScale->SetMarkerColor(kBlack);

  TF1* f1 = new TF1("line", "[0]", 200, 2100);
  f1->SetLineColor(kRed);
  gr_egmScale->Fit( f1, "QR" );

  gr_egmScale->Draw("Psame");

  ZGDrawTools::addLabels(c1, -1., "CMS Simulation");

  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/muEgmScaleSyst.eps", dir.c_str()));
  c1->SaveAs( Form("%s/muEgmScaleSyst.pdf", dir.c_str()));

  return 0;

}
