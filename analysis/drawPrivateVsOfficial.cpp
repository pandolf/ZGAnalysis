#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TLegend.h"


#include "../interface/ZGDrawTools.h"

void compareSingleMassPoint( TFile* file, int mass, const std::string& width );
void compareWithSelection( TTree* tree_private, TTree* tree_official, int mass, const std::string& width, const std::string& sel );


int main() {

  ZGDrawTools::setStyle();


  TFile* file = TFile::Open("EventYields_v0_eth74X/trees.root");

  compareSingleMassPoint( file, 750 , "0p014" );
  compareSingleMassPoint( file, 1500, "0p014" );
  compareSingleMassPoint( file, 2000, "0p014" );
  compareSingleMassPoint( file, 750 , "5p6" );
  compareSingleMassPoint( file, 1500, "5p6" );
  compareSingleMassPoint( file, 2000, "5p6" );


  return 0;

}


void compareSingleMassPoint( TFile* file, int mass, const std::string& width ) {

  TTree* tree_private = (TTree*)file->Get( Form("XZg_Spin0_ZToLL_W_%s_M_%d", width.c_str(), mass) );
  if( tree_private==0 ) {
    std::cout << "-> Didn't find private tree for mass point: " << mass << " (W=" << width << ")" << std::endl;
    return;
  }

  TTree* tree_official = (TTree*)file->Get( Form("GluGluSpin0ToZG_ZToLL_W%s_M%d", width.c_str(), mass) );
  if( tree_official==0 ) {
    std::cout << "-> Didn't find official tree for mass point: " << mass << " (W=" << width << ")" << std::endl;
    return;
  }


  compareWithSelection( tree_private, tree_official, mass, width, "ee" );
  compareWithSelection( tree_private, tree_official, mass, width, "mm" );

}


void compareWithSelection( TTree* tree_private, TTree* tree_official, int mass, const std::string& width, const std::string& sel ) {


  std::string tree_sel = (sel=="ee") ? "leptType==11" : "leptType==13";

  int nBins = 100;
  float xMin = 0.5*(float)mass;
  float xMax = 1.2*(float)mass;
  TH1D* h1_private = new TH1D( Form("mass_private_m%d_w%s_%s", mass, width.c_str(), sel.c_str()), "", nBins, xMin, xMax );
  h1_private->Sumw2();
  TH1D* h1_official = new TH1D( Form("mass_official_m%d_w%s_%s", mass, width.c_str(), sel.c_str()), "", nBins, xMin, xMax );
  h1_official->Sumw2();

  tree_private ->Project( h1_private->GetName() , "boss_mass", Form("weight*(%s)", tree_sel.c_str()) );
  tree_official->Project( h1_official->GetName(), "boss_mass", Form("weight*(%s)", tree_sel.c_str()) );

  float yMin = 0.005;
  float yHistoPriv = h1_private ->GetMaximum()/h1_private->Integral();
  float yHistoOffi = h1_official->GetMaximum()/h1_official->Integral();
  float yMaxHisto = (yHistoPriv>yHistoOffi) ? yHistoPriv : yHistoOffi;
  float yMax = 1.25*yMaxHisto;

  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();

  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0, yMax );
  h2_axes->SetXTitle("M(Z#gamma) [GeV]");
  h2_axes->SetYTitle("Normalized to Unity");
  h2_axes->Draw();


  h1_official->SetLineColor(kRed);
  h1_official->SetLineWidth(2);
  h1_official->SetMarkerColor(kRed);

  h1_private->SetLineWidth(2);
  
  h1_private ->DrawNormalized("same");
  h1_official->DrawNormalized("same");

  TString widthForLegend( width );
  widthForLegend.ReplaceAll("p", ".");
  std::string legendTitle = (sel=="ee") ? std::string(Form("M=%d, W=%s%%, ee#gamma", mass, widthForLegend.Data())) : std::string(Form("M=%d, W=%s%%, #mu#mu#gamma", mass, widthForLegend.Data()));

  TLegend* legend = new TLegend( 0.2, 0.68, 0.5, 0.92 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.035);
  legend->SetHeader( legendTitle.c_str() );
  legend->AddEntry( h1_private , "Private 74X" , "L" );
  legend->AddEntry( h1_official, "Official 76X", "L" );
  legend->Draw("same");

  gPad->RedrawAxis();

  ZGDrawTools::addLabels( c1, -1., "CMS Simulation" );

  c1->SaveAs( Form("privVsOff_m%d_w%s_%s.eps", mass, width.c_str(), sel.c_str()) );
  c1->SaveAs( Form("privVsOff_m%d_w%s_%s.pdf", mass, width.c_str(), sel.c_str()) );

  c1->Clear();
  c1->SetLogy();

  TH2D* h2_axes_log = new TH2D("axes_log", "", 10, xMin, xMax, 10, yMin, yMax );
  h2_axes_log->SetXTitle("M(Z#gamma) [GeV]");
  h2_axes_log->SetYTitle("Normalized to Unity");
  h2_axes_log->Draw();

  h1_private ->DrawNormalized("same");
  h1_official->DrawNormalized("same");

  legend->Draw("same");

  gPad->RedrawAxis();

  ZGDrawTools::addLabels( c1, -1., "CMS Simulation" );

  c1->SaveAs( Form("privVsOff_m%d_w%s_%s_log.eps", mass, width.c_str(), sel.c_str()) );
  c1->SaveAs( Form("privVsOff_m%d_w%s_%s_log.pdf", mass, width.c_str(), sel.c_str()) );

  delete c1;
  delete h2_axes;
  delete h2_axes_log;
  delete h1_private;
  delete h1_official;

}
