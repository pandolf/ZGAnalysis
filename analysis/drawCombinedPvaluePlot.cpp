#include <fstream>

#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"

#include "RooStats/RooStatsUtils.h"

#include "../interface/ZGDrawTools.h"
#include "../interface/ZGConfig.h"


void getPvalueGraphs( const std::string& pvaluesFile, TGraph* gr_pvalue );


int main( int argc, char* argv[] ) {


  if( argc<2 ) {
    std::cout << " USAGE: ./drawCombinedPvaluePlot [dir]" << std::endl;
    exit(1);
  }
  
  
  std::string dir(argv[1]);

  ZGDrawTools::setStyle();


  std::string pvaluesFile_comb  ( Form( "%s/pvalues_w0p014_comb_all.txt"   , dir.c_str() ));
  std::string pvaluesFile_only13( Form( "%s/pvalues_w0p014_comb_only13.txt", dir.c_str() ));
  std::string pvaluesFile_only8 ( Form( "%s/pvalues_w0p014_comb_only8.txt" , dir.c_str() ));


  std::string axisName;
  axisName = std::string("p_{0}");
  //axisName = std::string(Form("p-Value of the Background"));


  
  TGraph* gr_comb_obs = new TGraph(0);
  TGraph* gr_comb_exp = new TGraph(0);
  
  gr_comb_obs       ->SetName("comb_obs");
  gr_comb_exp       ->SetName("comb_exp");

  
  TGraph* gr_only13_obs = new TGraph(0);
  TGraph* gr_only13_exp = new TGraph(0);

  gr_only13_obs       ->SetName("only13_obs");
  gr_only13_exp       ->SetName("only13_exp");
  
  
  TGraph* gr_only8_obs = new TGraph(0);
  TGraph* gr_only8_exp = new TGraph(0);

  gr_only8_obs       ->SetName("only8_obs");
  gr_only8_exp       ->SetName("only8_exp");
  



  getPvalueGraphs( pvaluesFile_comb  , gr_comb_obs   );
  getPvalueGraphs( pvaluesFile_only13, gr_only13_obs );
  getPvalueGraphs( pvaluesFile_only8 , gr_only8_obs  );


  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();
  c1->SetLogy();
  c1->SetRightMargin(0.055);

  float xMin = 350.;
  float xMax = 2000.;
  float yMin = 0.001;

  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, yMin, 0.55 );
  h2_axes->SetYTitle( axisName.c_str() );
  //h2_axes->SetYTitle( "95\% CL UL on #sigma #times BR(A#rightarrowZ#gamma#rightarrowl^{+}l^{-}#gamma) [fb]");
  h2_axes->SetXTitle( "Resonance Mass [GeV]");
  h2_axes->GetYaxis()->SetTitleOffset(0.9);
  h2_axes->Draw();


  gr_only13_obs->SetLineStyle(2);
  gr_only8_obs ->SetLineStyle(2);

  gr_comb_obs  ->SetLineWidth(3);
  gr_only13_obs->SetLineWidth(3);
  gr_only8_obs ->SetLineWidth(3);

  gr_comb_obs  ->SetLineColor(kBlack);
  gr_only13_obs->SetLineColor(46);
  gr_only8_obs ->SetLineColor(8);

  gr_comb_exp  ->SetLineWidth(3);
  gr_only13_exp->SetLineWidth(3);
  gr_only8_exp ->SetLineWidth(3);

  gr_comb_exp  ->SetLineStyle(2);
  gr_only13_exp->SetLineStyle(2);
  gr_only8_exp ->SetLineStyle(2);

  gr_comb_exp  ->SetLineColor(kBlack);
  gr_only13_exp->SetLineColor(46);
  gr_only8_exp ->SetLineColor(8);


  TLegend* legend = new TLegend( 0.6, 0.2, 0.9, 0.45 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  legend->SetHeader("W = 0.014%");
  legend->AddEntry( gr_only8_obs , "8 TeV", "L" );
  legend->AddEntry( gr_only13_obs, "13 TeV", "L" );
  legend->AddEntry( gr_comb_obs, "Combination", "L" );
  legend->Draw("same");


  //TLegend* legend2 = new TLegend( 0.3, 0.75, 0.6, 0.9 );
  //legend2->SetFillColor(0);
  //legend2->SetTextSize(0.038);
  //legend2->SetTextFont(42);
  //legend2->AddEntry( gr_comb_exp, "Expected", "L" );
  //legend2->AddEntry( gr_comb_obs, "Observed", "L" );
  //legend2->Draw("same");

  //gr_only13_exp->Draw("L same");
  //gr_only8_exp ->Draw("L same");
  //gr_comb_exp  ->Draw("L same");



  for( unsigned i=1; i<5; i++ ) {

    float y = RooStats::SignificanceToPValue(i);

    TLine* line_sigma = new TLine( xMin, y, xMax, y );
    line_sigma->SetLineStyle(2);
  
    TLatex* text_sigma = new TLatex( xMax*1.02, y*0.9, Form("#color[%d]{%d#sigma}", kGray+2, i));
    text_sigma->SetTextSize(0.04);
    if( y>yMin ) {
      line_sigma->Draw("same");
      text_sigma->Draw("same");
    }

  }


  //float y_1sigma = 1.-0.6826894;
  //TLine* line_1sigma = new TLine( xMin, y_1sigma, xMax, y_1sigma );
  //line_1sigma->SetLineStyle(2);
  //line_1sigma->Draw("same");

  //TLatex* text_1sigma = new TLatex( xMax*1.01, y_1sigma*0.9, Form("#color[%d]{1#sigma}", kGray+2));
  //text_1sigma->SetTextSize(0.04);
  //text_1sigma->Draw("same");

  //TLine* line_2sigma = new TLine( xMin, 1.-0.954499, xMax, 1.-0.954499 );
  //line_2sigma->SetLineStyle(2);
  //line_2sigma->Draw("same");

  //TLine* line_3sigma = new TLine( xMin, 1.-0.997300, xMax, 1.-0.997300 );
  //line_3sigma->SetLineStyle(2);
  //line_3sigma->Draw("same");

  gr_only13_obs->Draw("L same");
  gr_only8_obs ->Draw("L same");
  gr_comb_obs  ->Draw("L same");



  ZGDrawTools::addLabels( c1, "CMS Preliminary, 19.7 fb^{-1} (8 TeV) + 2.7 fb^{-1} (13 TeV)");


  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/pvalue_comb_w0p014.eps" , dir.c_str() ));
  c1->SaveAs( Form("%s/pvalue_comb_w0p014.pdf" , dir.c_str() ));

  c1->SetLogx();
  h2_axes->GetXaxis()->SetMoreLogLabels();
  h2_axes->GetXaxis()->SetNoExponent();

  c1->SaveAs( Form("%s/pvalue_comb_w0p014_logx.eps" , dir.c_str() ));
  c1->SaveAs( Form("%s/pvalue_comb_w0p014_logx.pdf" , dir.c_str() ));



  return 0;

}



void getPvalueGraphs( const std::string& pvaluesFile, TGraph* gr_pvalue ) {

  ifstream ifs(pvaluesFile.c_str());
  std::cout << "-> Opened file: " << pvaluesFile << std::endl;
  int iPointExp = 0;
  int iPointObs = 0;
  float lastMass = -1;
  float lastObs = -1;

  while( ifs.good() ) {

    float m, pvalue;
    std::string s_m, s_pvalue;
    ifs >> s_m >> m >> s_pvalue >> pvalue;
    TString m_tstr(s_m);
    if( m_tstr.BeginsWith("#") ) continue;
    if( m==lastMass ) continue;
    //if( m==350. ) continue;


    std::cout << "m: " << m << " pvalue: " << pvalue << std::endl;


    if( pvalue>0. ) {
    //if( obs>0.001 && m!=610 && m!=620 && m!=730 && m!=920 && m!= 970 ) {
      gr_pvalue->SetPoint( iPointObs, m, pvalue );
      iPointObs++;
    }


    lastMass = m;
    //lastObs = obs;

  }

} 
