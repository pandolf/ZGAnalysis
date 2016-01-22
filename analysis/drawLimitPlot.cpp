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



int main( int argc, char* argv[] ) {


  if( argc==1 ) {
    std::cout << " USAGE: ./drawLimitPlot [cfg]" << std::endl;
    exit(1);
  }
  
  
  std::string configFileName(argv[1]);
  ZGConfig cfg(configFileName);


  ZGDrawTools::setStyle();

  std::string dir = cfg.getEventYieldDir() + "_fit_v0_default_shapes_combination_pcorr_lumi_2.3";
  std::string fullPath = "../../diphotons/Analysis/macros/" + dir;
  std::string limitsFile = fullPath + "/limits.txt";

  
  TGraph* gr_obs = new TGraph(0);
  TGraph* gr_exp = new TGraph(0);
  TGraphAsymmErrors* gr_exp_1sigma = new TGraphAsymmErrors(0);
  TGraphAsymmErrors* gr_exp_2sigma = new TGraphAsymmErrors(0);

  ifstream ifs(limitsFile.c_str());
  std::cout << "-> Opened file: " << limitsFile << std::endl;
  int iPoint = 0;

  while( ifs.good() ) {

    float m, obs, exp, exp_m1s, exp_m2s, exp_p1s, exp_p2s;
    std::string s_m, s_obs, s_exp, s_exp_m1s, s_exp_m2s, s_exp_p1s, s_exp_p2s;
    ifs >> s_m >> m >> s_obs >> obs >> s_exp >> exp >> s_exp_m1s >> exp_m1s >> s_exp_m2s >> exp_m2s >> s_exp_p1s >>  exp_p1s >> s_exp_p2s >> exp_p2s;
    if( m>=1000. ) continue;

    gr_obs       ->SetPoint( iPoint, m, obs );
    gr_exp       ->SetPoint( iPoint, m, exp );
    gr_exp_1sigma->SetPoint( iPoint, m, exp );
    gr_exp_2sigma->SetPoint( iPoint, m, exp );

    gr_exp_1sigma->SetPointError( iPoint, 0., 0., exp-exp_m1s, exp_p1s-exp );
    gr_exp_2sigma->SetPointError( iPoint, 0., 0., exp-exp_m2s, exp_p2s-exp );

    iPoint++;

  }

  gr_obs->SetLineWidth(2);

  gr_exp->SetLineWidth(2);
  gr_exp->SetLineStyle(2);

  gr_exp_1sigma->SetLineWidth(0);
  gr_exp_1sigma->SetFillColor(8);

  gr_exp_2sigma->SetLineWidth(0);
  gr_exp_2sigma->SetFillColor(219);

  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();

  TH2D* h2_axes = new TH2D("axes", "", 10, 350., 950., 10, 0., 55. );
  h2_axes->SetYTitle( "95\% CL Upper Limit [fb]");
  h2_axes->SetXTitle( "Resonance Mass [GeV]");
  h2_axes->Draw();


  gr_exp_2sigma->Draw("E3 same");
  gr_exp_1sigma->Draw("E3 same");
  gr_exp       ->Draw("L  same");
  //gr_obs       ->Draw("L  same");


  gr_exp_1sigma->SetLineWidth(2);
  gr_exp_1sigma->SetLineStyle(2);
  gr_exp_2sigma->SetLineWidth(2);
  gr_exp_2sigma->SetLineStyle(2);

  TLegend* legend = new TLegend( 0.55, 0.7, 0.9, 0.9 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  //legend->AddEntry( gr_exp, "Expected", "L" );
  legend->AddEntry( gr_exp_1sigma, "Expected #pm 1#sigma", "LF" );
  legend->AddEntry( gr_exp_2sigma, "Expected #pm 2#sigma", "LF" );
  //legend->AddEntry( gr_obs, "Observed", "L" );
  legend->Draw("same");


  ZGDrawTools::addLabels( c1, 2.3, "CMS Simulation");


  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/limit.eps", cfg.getEventYieldDir().c_str()) );
  c1->SaveAs( Form("%s/limit.pdf", cfg.getEventYieldDir().c_str()) );

  return 0;

}
