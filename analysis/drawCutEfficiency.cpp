#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

#include "../interface/ZGDrawTools.h"
#include "../interface/ZGConfig.h"


void drawEfficiency( const std::string& outdir, TTree* tree, const std::string& saveName, const std::string& cut, const std::string& labelText );


int main( int argc, char* argv[] ) {


  if( argc==1 ) {
    std::cout << " USAGE: ./drawCutEfficiency [cfg] [lumi/shape]" << std::endl;
    exit(1);
  }
  
  std::string configFileName(argv[1]);
  ZGConfig cfg(configFileName);

  ZGDrawTools::setStyle();

  std::string theDir = cfg.getEventYieldDir();

  TFile* file = TFile::Open(Form("%s/trees.root", theDir.c_str()));
  TTree* tree = (TTree*)file->Get("zg");

  drawEfficiency( theDir, tree, "ptCut", "gamma_pt/boss_mass > 40./150.", "p_{T}^{#gamma} > #frac{40}{150} M(Z#gamma)" );
  drawEfficiency( theDir, tree, "etaCut", "abs(gamma_eta)<1.44", "|#eta^{#gamma}| < 1.44" );

  return 0;

}


void drawEfficiency( const std::string& outdir, TTree* tree, const std::string& saveName, const std::string& cut, const std::string& labelText ) {


  int nBins = 8;
  Double_t bins[nBins+1];
  bins[0] = 300.;
  bins[1] = 350.;
  bins[2] = 400.;
  bins[3] = 450.;
  bins[4] = 500.;
  bins[5] = 600.;
  bins[6] = 700.;
  bins[7] = 800.;
  bins[8] = 950.;


  TH1D* h1_denom = new TH1D("denom", "", nBins, bins );
  h1_denom->Sumw2();

  TH1D* h1_num = new TH1D("num", "", nBins, bins );
  h1_num->Sumw2();

  tree->Project( "denom", "boss_mass", "weight" );
  tree->Project( "num", "boss_mass", Form("weight*(%s)", cut.c_str()) );

  h1_num->Divide( h1_denom );

  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();

  TH2D* h2_axes = new TH2D("axes", "", 10, bins[0], bins[nBins], 10, 0., 1.0001 );
  h2_axes->SetXTitle("M(Z#gamma) [GeV]" );
  h2_axes->SetYTitle("Efficiency");
  h2_axes->Draw();

  TPaveText* label = new TPaveText( 0.2, 0.25, 0.55, 0.3, "brNDC" );
  label->SetTextSize(0.038 );
  label->SetTextFont(42);
  label->SetFillColor(0);
  label->AddText( labelText.c_str() );
  label->Draw("same");


  TF1* line = new TF1( "f1", "[0] + [1]*x", bins[0], bins[nBins] );
  h1_num->Fit( line, "QR0" );
  TH1D* band = ZGDrawTools::getBand( line );
  band->Draw("C E3 same" );
  line->SetLineColor(46);
  line->Draw("same");

  h1_num->SetMarkerStyle(20);
  h1_num->SetMarkerSize(1.3);
  h1_num->Draw("p same");


  gPad->RedrawAxis();

  ZGDrawTools::addLabels(c1, -1, "CMS Simulation");

  c1->SaveAs( Form("%s/eff_%s.eps", outdir.c_str(), saveName.c_str() ) );
  c1->SaveAs( Form("%s/eff_%s.pdf", outdir.c_str(), saveName.c_str() ) );

  delete c1;
  delete h2_axes;
  delete h1_num;
  delete h1_denom;
  delete band;
  delete line;

}
