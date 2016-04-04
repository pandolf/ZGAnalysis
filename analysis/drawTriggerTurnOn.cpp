#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TLegend.h"

#include "../interface/ZGDrawTools.h"

TEfficiency* getTurnOn( TTree* tree, const std::string& hlt, const std::string& selection="" );


float xMin = 200.;
float xMax = 1000.;



int main() {

  ZGDrawTools::setStyle();

  TFile* file = TFile::Open("EventYields_v0_eth74X_noHLT/trees.root");
  TTree* tree_zg = (TTree*)file->Get("zg");

  std::vector<std::string> selections;
  selections.push_back("ee");
  selections.push_back("mm");

  std::vector<TEfficiency*> v_eff;

  for( unsigned i=0; i<selections.size(); ++i ) {

    std::string trig1, trig2;
    if( selections[i]=="ee" ) {
      trig1 = "HLT_DoubleEle";
      trig2 = "HLT_DoubleEle33";
    } else if( selections[i]=="mm" ) {
      trig1 = "HLT_DoubleMu";
      trig2 = "HLT_SingleMu";
    }
    TEfficiency* eff_trig1 = getTurnOn( tree_zg, trig1, selections[i] );
    TEfficiency* eff_trig2 = getTurnOn( tree_zg, trig2, selections[i] );
    TEfficiency* eff_all   = getTurnOn( tree_zg, trig1 + " || " + trig2, selections[i] );
    eff_all->SetName(Form("eff_all_%s", selections[i].c_str()));

    TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
    c1->cd();

    eff_trig1->SetMarkerSize(1.1);
    eff_trig1->SetMarkerColor(kBlack);
    eff_trig1->SetMarkerStyle(20);

    eff_trig2->SetMarkerSize(1.1);
    eff_trig2->SetMarkerColor(kRed);
    eff_trig2->SetLineColor(kRed);
    eff_trig2->SetMarkerStyle(20);

    eff_all->SetMarkerSize(1.1);
    eff_all->SetMarkerColor(kGreen);
    eff_all->SetLineColor(kGreen);
    eff_all->SetMarkerStyle(20);


    TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., 1.1 );
    h2_axes->SetXTitle( "M(Z#gamma) [GeV]");
    h2_axes->SetYTitle( "Efficiency" );
    h2_axes->Draw();

    TLine* lineone = new TLine( xMin, 1., xMax, 1. );
    lineone->Draw("same");

    eff_trig1->Draw("same");
    eff_trig2->Draw("same");
    eff_all->Draw("same");

    TLegend* legend = new TLegend( 0.55, 0.2, 0.9, 0.46 );
    legend->SetFillColor(0);
    legend->SetTextSize(0.035);
    if( selections[i]=="ee" ) {
      legend->SetHeader( "Electron Channel" );
      eff_all->SetTitle( "Electron Channel" );
    } else if( selections[i]=="mm" ) {
      legend->SetHeader( "Muon Channel" );
      eff_all->SetTitle( "Muon Channel" );
    }
    legend->AddEntry( eff_trig1, trig1.c_str(), "P" );
    legend->AddEntry( eff_trig2, trig2.c_str(), "P" );
    legend->AddEntry( eff_all, "Both in OR", "P" );
    legend->Draw("same");

    ZGDrawTools::addLabels( c1, -1, "CMS Simulation" );

    gPad->RedrawAxis();

    c1->SaveAs( Form("turnOn_%s.eps", selections[i].c_str()) );
    c1->SaveAs( Form("turnOn_%s.pdf", selections[i].c_str()) );

    delete c1;
    delete h2_axes;

    v_eff.push_back( eff_all );

  }

  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();

  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., 1.1 );
  h2_axes->SetXTitle( "M(Z#gamma) [GeV]");
  h2_axes->SetYTitle( "Efficiency" );
  h2_axes->Draw();

  TLine* lineone = new TLine( xMin, 1., xMax, 1. );
  lineone->Draw("same");

  TLegend* legend = new TLegend( 0.2, 0.2, 0.5, 0.4 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.035);

  for( unsigned i=0; i<v_eff.size(); ++i ) {
    v_eff[i]->SetMarkerColor(i+1);
    v_eff[i]->SetLineColor(i+1);
    v_eff[i]->Draw("p same");
    legend->AddEntry( v_eff[i], v_eff[i]->GetTitle(), "P" );
  }

  legend->Draw("same");

  ZGDrawTools::addLabels( c1, -1, "CMS Simulation" );

  gPad->RedrawAxis();

  c1->SaveAs( "turnOn_ee_vs_mm.eps" );
  c1->SaveAs( "turnOn_ee_vs_mm.pdf" );


  return 0;

}


TEfficiency* getTurnOn( TTree* tree, const std::string& hlt, const std::string& selection ) {

  int nBins = 25;

  TH1D* h1_num   = new TH1D(Form("num_%s_%s"  , hlt.c_str(), selection.c_str()), "", nBins, xMin, xMax );
  h1_num->Sumw2();
  TH1D* h1_denom = new TH1D(Form("denom_%s_%s", hlt.c_str(), selection.c_str()), "", nBins, xMin, xMax );
  h1_denom->Sumw2();


  std::string selTree = "";
  if( selection=="ee" )
    selTree="leptType==11";
  else if( selection=="mm" )
    selTree="leptType==13";
  else if( selection=="all" )
    selTree="boss_mass>0.";

  tree->Project( h1_num  ->GetName(), "boss_mass", Form("(%s) && (%s)", hlt.c_str(), selTree.c_str()) );
  tree->Project( h1_denom->GetName(), "boss_mass", selTree.c_str() );

  TEfficiency* turnon = new TEfficiency( *h1_num, *h1_denom );
  turnon->SetName( Form("eff_%s_%s", hlt.c_str(), selection.c_str()) );

  return turnon;

}
