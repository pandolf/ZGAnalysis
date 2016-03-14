#include <cmath>

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TF1.h"
#include "TLegend.h"


#include "../interface/ZGDrawTools.h"

void compareSingleMassPoint( TFile* file, int mass, const std::string& width, TGraphErrors* eff_ee, TGraphErrors* eff_mm );
void compareWithSelection( TTree* tree_private, TTree* tree_official, int mass, const std::string& width, const std::string& sel, TGraphErrors* eff );
void compareEffToPriv( TGraphErrors* eff, TF1* f1_eff_priv, const std::string& width, const std::string& sel );


int main() {

  ZGDrawTools::setStyle();


  system("mkdir -p privateVsOfficial");

  TFile* file = TFile::Open("EventYields_v0_eth74X/trees.root");

  TGraphErrors* eff_w0p014_ee = new TGraphErrors(0);
  TGraphErrors* eff_w0p014_mm = new TGraphErrors(0);
  TGraphErrors* eff_w5p6_ee   = new TGraphErrors(0);
  TGraphErrors* eff_w5p6_mm   = new TGraphErrors(0);

  eff_w0p014_ee->SetName("eff_w0p014_ee");
  eff_w0p014_mm->SetName("eff_w0p014_mm");
  eff_w5p6_ee  ->SetName("eff_w5p6_ee");
  eff_w5p6_mm  ->SetName("eff_w5p6_mm");

  compareSingleMassPoint( file, 750 , "0p014", eff_w0p014_ee, eff_w0p014_mm );
  compareSingleMassPoint( file, 1500, "0p014", eff_w0p014_ee, eff_w0p014_mm );
  compareSingleMassPoint( file, 2000, "0p014", eff_w0p014_ee, eff_w0p014_mm );
  compareSingleMassPoint( file, 750 , "5p6"  , eff_w5p6_ee  , eff_w5p6_mm   );
  compareSingleMassPoint( file, 1500, "5p6"  , eff_w5p6_ee  , eff_w5p6_mm   );
  compareSingleMassPoint( file, 2000, "5p6"  , eff_w5p6_ee  , eff_w5p6_mm   );


  // compare to private efficiency:
  TFile* effFile_w0p014_private = TFile::Open("EventYields_v0_eth74X/signalEfficiency_w0p014.root");
  TFile* effFile_w5p6_private   = TFile::Open("EventYields_v0_eth74X/signalEfficiency_w5p6.root");

  TF1* f1_priv_w0p014_ee = (TF1*)effFile_w0p014_private->Get("f1_eff_ee");
  TF1* f1_priv_w0p014_mm = (TF1*)effFile_w0p014_private->Get("f1_eff_mm");
  TF1* f1_priv_w5p6_ee   = (TF1*)effFile_w5p6_private  ->Get("f1_eff_ee");
  TF1* f1_priv_w5p6_mm   = (TF1*)effFile_w5p6_private  ->Get("f1_eff_mm");

  compareEffToPriv( eff_w0p014_ee, f1_priv_w0p014_ee, "0p014", "ee" );
  compareEffToPriv( eff_w0p014_mm, f1_priv_w0p014_mm, "0p014", "mm" );
  compareEffToPriv( eff_w5p6_ee  , f1_priv_w5p6_ee  , "5p6"  , "ee" );
  compareEffToPriv( eff_w5p6_mm  , f1_priv_w5p6_mm  , "5p6"  , "mm" );

  TFile* effFile = TFile::Open( "privateVsOfficial/eff_official.root", "recreate" );
  effFile->cd();
  eff_w0p014_ee->Write();
  eff_w0p014_mm->Write();
  eff_w5p6_ee  ->Write();
  eff_w5p6_mm  ->Write();
  effFile->Close();

  std::cout << "-> Saved official efficiencies in: " << effFile->GetName() << std::endl;

  return 0;

}


void compareSingleMassPoint( TFile* file, int mass, const std::string& width, TGraphErrors* eff_ee, TGraphErrors* eff_mm ) {

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


  compareWithSelection( tree_private, tree_official, mass, width, "ee", eff_ee );
  compareWithSelection( tree_private, tree_official, mass, width, "mm", eff_mm );

}


void compareWithSelection( TTree* tree_private, TTree* tree_official, int mass, const std::string& width, const std::string& sel, TGraphErrors* eff ) {


  std::string tree_sel = (sel=="ee") ? "leptType==11" : "leptType==13";

  // firs off efficiency (for official only)
  Double_t denom = 50000./3.;
  Double_t num = (Double_t)tree_official->GetEntries(tree_sel.c_str());
  Double_t thisEff = num/denom;
  Double_t thisEff_err = sqrt( thisEff*(1.-thisEff)/denom );
  int iPoint = eff->GetN();
  eff->SetPoint     (iPoint, mass, thisEff    );
  eff->SetPointError(iPoint,   0., thisEff_err);

  int nBins = 100;
  float xMin = 0.65*(float)mass;
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

  c1->SaveAs( Form("privateVsOfficial/privVsOff_m%d_w%s_%s.eps", mass, width.c_str(), sel.c_str()) );
  c1->SaveAs( Form("privateVsOfficial/privVsOff_m%d_w%s_%s.pdf", mass, width.c_str(), sel.c_str()) );

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

  c1->SaveAs( Form("privateVsOfficial/privVsOff_m%d_w%s_%s_log.eps", mass, width.c_str(), sel.c_str()) );
  c1->SaveAs( Form("privateVsOfficial/privVsOff_m%d_w%s_%s_log.pdf", mass, width.c_str(), sel.c_str()) );

  delete c1;
  delete h2_axes;
  delete h2_axes_log;
  delete h1_private;
  delete h1_official;

  

}


void compareEffToPriv( TGraphErrors* eff, TF1* f1_eff_priv, const std::string& width, const std::string& sel ) {


  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  TH2D* h2_axes = new TH2D("axes", "", 10, 300., 2100., 10, 0., 1.);
  h2_axes->SetXTitle( "M(Z#gamma) [GeV]");
  h2_axes->SetYTitle( "Efficiency" );
  h2_axes->Draw("same");

  f1_eff_priv->SetLineColor(46);
  f1_eff_priv->Draw("same");
  
  eff->SetMarkerStyle(20);
  eff->SetMarkerSize(1.5);
  eff->Draw("p same");

  TString widthForLegend( width );
  widthForLegend.ReplaceAll("p", ".");
  std::string legendTitle = (sel=="ee") ? std::string(Form("W=%s%%, ee#gamma", widthForLegend.Data())) : std::string(Form("W=%s%%, #mu#mu#gamma", widthForLegend.Data()));

  TLegend* legend = new TLegend( 0.2, 0.68, 0.5, 0.92 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.035);
  legend->SetHeader( legendTitle.c_str() );
  legend->AddEntry( f1_eff_priv , "From Private Samples" , "L" );
  legend->AddEntry( eff, "From Official 76X", "P" );
  legend->Draw("same");

  gPad->RedrawAxis();

  ZGDrawTools::addLabels( c1, -1., "CMS Simulation" );

  c1->SaveAs( Form("privateVsOfficial/eff_privVsOff_w%s_%s.eps", width.c_str(), sel.c_str()) );
  c1->SaveAs( Form("privateVsOfficial/eff_privVsOff_w%s_%s.pdf", width.c_str(), sel.c_str()) );
  
  delete h2_axes;
  delete c1;

}
