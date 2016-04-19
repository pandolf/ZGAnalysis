#include <fstream>

#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphAsymmErrors.h"

#include "../interface/ZGDrawTools.h"
#include "../interface/ZGConfig.h"


void getLimitGraphs( const std::string& limitsFile, TGraph* gr_obs, TGraph* gr_exp, TGraphAsymmErrors* gr_exp_1sigma, TGraphAsymmErrors* gr_exp_2sigma );


int main( int argc, char* argv[] ) {


  if( argc<2 ) {
    std::cout << " USAGE: ./drawCombinedLimitPlot [dir]" << std::endl;
    exit(1);
  }
  
  
  std::string dir(argv[1]);

  ZGDrawTools::setStyle();


  std::string limitsFile_comb  ( Form( "%s/limits_w0p014_comb_all.txt"   , dir.c_str() ));
  std::string limitsFile_only13( Form( "%s/limits_w0p014_comb_only13.txt", dir.c_str() ));
  std::string limitsFile_only8 ( Form( "%s/limits_w0p014_comb_only8.txt" , dir.c_str() ));


  std::string axisName;
  axisName = std::string(Form("95\%% CL UL on #sigma #times BR(A#rightarrowZ#gamma#rightarrow l^{+}l^{-}#gamma) [fb]"));


  
  TGraph* gr_comb_obs = new TGraph(0);
  TGraph* gr_comb_exp = new TGraph(0);
  TGraphAsymmErrors* gr_comb_exp_1sigma = new TGraphAsymmErrors(0);
  TGraphAsymmErrors* gr_comb_exp_2sigma = new TGraphAsymmErrors(0);
  
  TGraph* gr_only13_obs = new TGraph(0);
  TGraph* gr_only13_exp = new TGraph(0);
  TGraphAsymmErrors* gr_only13_exp_1sigma = new TGraphAsymmErrors(0);
  TGraphAsymmErrors* gr_only13_exp_2sigma = new TGraphAsymmErrors(0);
  
  TGraph* gr_only8_obs = new TGraph(0);
  TGraph* gr_only8_exp = new TGraph(0);
  TGraphAsymmErrors* gr_only8_exp_1sigma = new TGraphAsymmErrors(0);
  TGraphAsymmErrors* gr_only8_exp_2sigma = new TGraphAsymmErrors(0);

  getLimitGraphs( limitsFile_comb, gr_comb_obs, gr_comb_exp, gr_comb_exp_1sigma, gr_comb_exp_2sigma );
  getLimitGraphs( limitsFile_only13, gr_only13_obs, gr_only13_exp, gr_only13_exp_1sigma, gr_only13_exp_2sigma );
  getLimitGraphs( limitsFile_only8, gr_only8_obs, gr_only8_exp, gr_only8_exp_1sigma, gr_only8_exp_2sigma );


  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();

  TH2D* h2_axes = new TH2D("axes", "", 10, 350., 1200., 10, 0., 9. );
  h2_axes->SetYTitle( axisName.c_str() );
  //h2_axes->SetYTitle( "95\% CL UL on #sigma #times BR(A#rightarrowZ#gamma#rightarrowl^{+}l^{-}#gamma) [fb]");
  h2_axes->SetXTitle( "Resonance Mass [GeV]");
  h2_axes->Draw();


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


  TLegend* legend = new TLegend( 0.6, 0.65, 0.9, 0.9 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  legend->SetHeader("W = 0.014%");
  legend->AddEntry( gr_only8_obs , "8 TeV", "L" );
  legend->AddEntry( gr_only13_obs, "13 TeV", "L" );
  legend->AddEntry( gr_comb_obs, "Combination", "L" );
  legend->Draw("same");


  TLegend* legend2 = new TLegend( 0.3, 0.75, 0.6, 0.9 );
  legend2->SetFillColor(0);
  legend2->SetTextSize(0.038);
  legend2->SetTextFont(42);
  legend2->AddEntry( gr_comb_exp, "Expected", "L" );
  legend2->AddEntry( gr_comb_obs, "Observed", "L" );
  legend2->Draw("same");

  gr_only13_exp->Draw("L same");
  gr_only8_exp ->Draw("L same");
  gr_comb_exp  ->Draw("L same");

  gr_only13_obs->Draw("L same");
  gr_only8_obs ->Draw("L same");
  gr_comb_obs  ->Draw("L same");



  ZGDrawTools::addLabels( c1, "CMS Preliminary, 19.7 fb^{-1} (8 TeV) + 2.7 fb^{-1} (13 TeV)");


  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/limit_comb_w0p014.eps" , dir.c_str() ));
  c1->SaveAs( Form("%s/limit_comb_w0p014.pdf" , dir.c_str() ));

  return 0;

}



void getLimitGraphs( const std::string& limitsFile, TGraph* gr_obs, TGraph* gr_exp, TGraphAsymmErrors* gr_exp_1sigma, TGraphAsymmErrors* gr_exp_2sigma ) {

  ifstream ifs(limitsFile.c_str());
  std::cout << "-> Opened file: " << limitsFile << std::endl;
  int iPointExp = 0;
  int iPointObs = 0;
  float lastMass = -1;
  float lastObs = -1;

  while( ifs.good() ) {

    float m, obs, exp, exp_m1s, exp_m2s, exp_p1s, exp_p2s;
    std::string s_m, s_obs, s_exp, s_exp_m1s, s_exp_m2s, s_exp_p1s, s_exp_p2s;
    ifs >> s_m >> m >> s_obs >> obs >> s_exp >> exp >> s_exp_m1s >> exp_m1s >> s_exp_m2s >> exp_m2s >> s_exp_p1s >>  exp_p1s >> s_exp_p2s >> exp_p2s;
    TString m_tstr(s_m);
    if( m_tstr.BeginsWith("#") ) continue;
    if( m==lastMass ) continue;
    //if( m==350. ) continue;


    //std::cout << "m: " << m << " obs: " << obs << " exp: " << exp << " exp_m1s: " << exp_m1s << " exp_m2s: " << exp_m2s << " exp_p1s: " << exp_p1s << " exp_p2s: " << exp_p2s << std::endl;



    if( obs>0.001 && m!=610 && m!=620 && m!=730 && m!=920 && m!= 970 ) {
      gr_obs       ->SetPoint( iPointObs, m, obs );
      iPointObs++;
    }

    //bool okForExp = ((m==348.) || (m==2000.) || (m==748.)
    //              || (m>=400. && m<900.  && (int(m) % 50 == 0) && m!=750.) 
    //              //|| (m>=700. && m<1100. && (int(m) % 100 == 0))
    //              || (m>=900. && (int(m-100) % 200 == 0))  );
    //if( okForExp ) {
    gr_exp       ->SetPoint( iPointExp, m, exp );
    gr_exp_1sigma->SetPoint( iPointExp, m, exp );
    gr_exp_2sigma->SetPoint( iPointExp, m, exp );
    gr_exp_1sigma->SetPointError( iPointExp, 0., 0., exp-exp_m1s, exp_p1s-exp );
    gr_exp_2sigma->SetPointError( iPointExp, 0., 0., exp-exp_m2s, exp_p2s-exp );
    iPointExp++;

    lastMass = m;
    lastObs = obs;

  }

} 
