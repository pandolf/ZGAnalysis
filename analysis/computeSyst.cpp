#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TProfile.h"

#include "../interface/ZGDrawTools.h"
#include "../interface/ZGConfig.h"



void computeSyst( const std::string& outputdir, TTree* tree, const std::string& weightvar, const std::string& systName );



int main( int argc, char* argv[] ) {

  if( argc==1 ) {
    std::cout << " USAGE: ./computeSignalSyst [cfg]" << std::endl;
    exit(1);
  }
  
  std::string configFileName(argv[1]);
  ZGConfig cfg(configFileName);

  std::string treesFile = cfg.getEventYieldDir() + "/trees.root";
  TFile* file = TFile::Open(treesFile.c_str());
  TTree* tree = (TTree*)file->Get("zg");

  ZGDrawTools::setStyle();

  std::string outputdir = cfg.getEventYieldDir();

  std::cout << "-> Starting PDF..." << std::endl;
  computeSyst( outputdir, tree, "weight_pdf", "PDF" );
  std::cout << "-> Starting Scale..." << std::endl;
  computeSyst( outputdir, tree, "weight_scale", "Scale" );
  std::cout << "DONE." << std::endl;


  return 0;

}



void computeSyst( const std::string& outputdir, TTree* tree, const std::string& weightvar, const std::string& systName ) {

  int nBins = 7;
  Double_t bins[nBins+1];
  bins[0] = 350.;
  bins[1] = 400.;
  bins[2] = 450.;
  bins[3] = 500.;
  bins[4] = 600.;
  bins[5] = 700.;
  bins[6] = 800.;
  bins[7] = 950.;

  TH1D* h1_ref = new TH1D("ref", "", nBins, bins );
  h1_ref->Sumw2();

  TH1D* h1_syst = new TH1D("syst", "", nBins, bins );
  h1_syst->Sumw2();


  tree->Project( "ref" , "boss_mass", "weight" );
  tree->Project( "syst", "boss_mass", weightvar.c_str() );


  h1_syst->Divide( h1_ref );
  for( int iBin=1; iBin<h1_syst->GetXaxis()->GetNbins()+1; ++iBin ) {
    h1_syst->SetBinContent(iBin, h1_syst->GetBinContent(iBin)-1. );
    h1_syst->SetBinError(iBin, h1_ref->GetBinError(iBin)/h1_ref->GetBinContent(iBin)*h1_syst->GetBinContent(iBin) );
  }

  h1_syst->Scale(100.); // in percent

  TCanvas* c1 = new TCanvas( "c2", "c2", 600, 600 );
  c1->cd();

  TF1* line = new TF1( "f1", "[0] + [1]*x", bins[0], bins[nBins] );
  h1_syst->Fit( line, "QR0" );
  TH1D* band = ZGDrawTools::getBand( line );


  TH2D* h2_axes = new TH2D("axes", "", 10, bins[0], bins[nBins], 10, 0., 20. );
  h2_axes->SetXTitle("M(Z#gamma) [GeV]");
  h2_axes->SetYTitle( Form("%s Systematic [%%]", systName.c_str() ) );
  h2_axes->Draw();

  h1_syst->SetMarkerStyle( 20 );
  h1_syst->SetMarkerColor( 46 );
  h1_syst->SetLineColor( 46 );
  h1_syst->SetMarkerSize( 1.3 );


  band->Draw("C E3 same" );
  line->Draw("same");
  h1_syst->Draw("p same");

  ZGDrawTools::addLabels( c1, -1., "CMS Simulation" );

  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/syst%s.eps", outputdir.c_str(), systName.c_str()) );
  c1->SaveAs( Form("%s/syst%s.pdf", outputdir.c_str(), systName.c_str()) );

  delete h2_axes;
  delete c1;
  delete line;
  delete h1_ref;
  delete h1_syst;
  delete band;

}
