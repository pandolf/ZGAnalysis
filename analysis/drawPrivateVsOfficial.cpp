#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLegend.h"


#include "../interface/ZGDrawTools.h"


void drawVar( TTree* tree_official, TTree* tree_private, const std::string& savename, const std::string& name, int nBins, float xMin, float xMax, const std::string& selection="" );

int main() {

  ZGDrawTools::setStyle();

  TFile* file = TFile::Open("EventYields_v0_2016v5/trees.root");

  std::vector<std::string> widths;
  //widths.push_back( "0p014" );
  widths.push_back( "5p6" );

  for( unsigned i=0; i<widths.size(); ++i ) {

    std::string name_official = "GluGluSpin0ToZG_ZToLL_W" + widths[i] + "_M500";
    std::string name_private  = "XZg_Spin0ToZG_ZToLL_W_" + widths[i] + "_M_500";

    TTree* tree_official = (TTree*)file->Get(name_official.c_str());
    TTree* tree_private  = (TTree*)file->Get(name_private.c_str());

    drawVar( tree_official, tree_private, "mZg_ee_" + widths[i], "boss_mass", 60, 250., 750., "leptType==11" );
    drawVar( tree_official, tree_private, "mZg_mm_" + widths[i], "boss_mass", 60, 250., 750., "leptType==13" );
  //drawVar( tree_official, tree_private, "gamma_pt","gamma_pt", 25, 40., 300.);
  //drawVar( tree_official, tree_private, "mZee", "z_mass", 80, 50., 130., "leptType==11");
  //drawVar( tree_official, tree_private, "mZmm", "z_mass", 80, 50., 130., "leptType==13");

  }

  return 0;

}

void drawVar( TTree* tree_official, TTree* tree_private, const std::string& savename, const std::string& name, int nBins, float xMin, float xMax, const std::string& selection ) {

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();
 

  float binWidth = (xMax-xMin)/((float)nBins);

  TH1D* h1_official = new TH1D( "h1_official", "", nBins, xMin, xMax );
  TH1D* h1_private = new TH1D( "h1_private", "", nBins, xMin, xMax );

  h1_official->Sumw2();
  h1_private ->Sumw2();

  h1_private->SetXTitle( name.c_str() );
  //h1_private->SetXTitle( "M(Z#gamma) [GeV]");
  h1_private->SetYTitle( "Events" );
  //h1_private->SetYTitle( "Normalized to Unity" );
  //h1_private->SetYTitle( Form("Events / (%.0f GeV)", binWidth) );

  h1_official->SetLineWidth(2);
  h1_private->SetLineWidth(2);

  if( selection=="" ) {
    tree_official->Project("h1_official", name.c_str(), "weight" );
    tree_private ->Project("h1_private" , name.c_str(), "weight" );
  } else {
    tree_official->Project("h1_official", name.c_str(), Form("weight*(%s)", selection.c_str()) );
    tree_private ->Project("h1_private" , name.c_str(), Form("weight*(%s)", selection.c_str()) );
  }
  h1_private->SetLineColor(kRed);
  h1_private->Draw("same");
  h1_official->Draw("same");
  //h1_private->DrawNormalized("same");
  //h1_official->DrawNormalized("same");

  ZGDrawTools::addLabels( c1, -1., "CMS Simulation" );

  TLegend* legend = new TLegend( 0.6, 0.75, 0.9, 0.9 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  legend->AddEntry( h1_official, "Official Sample", "L" );
  legend->AddEntry( h1_private, "Private Sample", "L" );
  legend->Draw("same");

  gPad->RedrawAxis();

  c1->SaveAs( Form("%s_privVsOff.eps", savename.c_str()) );
  c1->SaveAs( Form("%s_privVsOff.pdf", savename.c_str()) );

  c1->SetLogy();

  c1->SaveAs( Form("%s_privVsOff_log.eps", savename.c_str()) );
  c1->SaveAs( Form("%s_privVsOff_log.pdf", savename.c_str()) );
  
  delete c1;
  delete h1_official;
  delete h1_private;

}
