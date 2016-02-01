#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TLegend.h"

#include "../interface/ZGDrawTools.h"




int main() {

  ZGDrawTools::setStyle();

  TFile* file = TFile::Open( "EventYields_v0/trees.root" );
  TTree* tree = (TTree*)file->Get("zg");

  float xMin = 200.;
  float xMax = 2200.;
  TH1D* h1_mass_pos = new TH1D("mass_pos", "", 40, xMin, xMax );
  TH1D* h1_mass_neg = new TH1D("mass_neg", "", 40, xMin, xMax );
  h1_mass_pos->Sumw2();
  h1_mass_neg->Sumw2();
  
  h1_mass_pos->SetLineWidth(2);
  h1_mass_neg->SetLineWidth(2);

  h1_mass_pos->SetLineColor(46);
  h1_mass_neg->SetLineColor(38);
  
  tree->Project( "mass_pos", "boss_mass", "weight>0." );
  tree->Project( "mass_neg", "boss_mass", "weight<0." );

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  float yMax = h1_mass_pos->GetMaximum()/h1_mass_pos->Integral()*1.2;
  
  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, 1200., 10, 0., yMax );
  h2_axes->SetXTitle( "M(Z#gamma) [GeV]" );
  h2_axes->SetYTitle( "Normalized to Unity" );
  h2_axes->Draw("");

  h1_mass_pos->DrawNormalized("same");
  h1_mass_neg->DrawNormalized("same");

  TLegend* legend = new TLegend( 0.6, 0.7, 0.9, 0.9 );
  legend->SetFillColor(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.038);
  legend->AddEntry( h1_mass_pos, "w > 0", "L" );
  legend->AddEntry( h1_mass_neg, "w < 0", "L" );
  legend->Draw("same");

  ZGDrawTools::addLabels( c1, -1, "CMS Simulation" );

  gPad->RedrawAxis();

  c1->SaveAs("massWeightCompare.eps");
  c1->SaveAs("massWeightCompare.pdf");

  c1->Clear();
  c1->SetLogy();

  TH2D* h2_axes_log = new TH2D("axes_log", "", 10, xMin, xMax, 10, 0.0001, 2.*yMax );
  h2_axes_log->SetXTitle( "M(Z#gamma) [GeV]" );
  h2_axes_log->SetYTitle( "Normalized to Unity" );
  h2_axes_log->Draw("");

  h1_mass_pos->Rebin(2);
  h1_mass_neg->Rebin(2);

  h1_mass_pos->DrawNormalized("same");
  h1_mass_neg->DrawNormalized("same");

  legend->Draw("same");

  ZGDrawTools::addLabels( c1, -1, "CMS Simulation" );

  gPad->RedrawAxis();
  
  c1->SaveAs("massWeightCompare_log.eps");
  c1->SaveAs("massWeightCompare_log.pdf");

  return 0;

}
