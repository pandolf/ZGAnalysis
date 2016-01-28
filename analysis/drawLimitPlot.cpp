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
    std::cout << " USAGE: ./drawLimitPlot [cfg] [eff=0.4]" << std::endl;
    exit(1);
  }
  
  
  std::string configFileName(argv[1]);
  ZGConfig cfg(configFileName);

  int eff = 40;
  if( argc>2 ) {
    std::string eff_str(argv[2]);
    float eff_float = atof(eff_str.c_str());
    if( eff_float<1. ) eff_float *= 100.;
    eff = eff_float;
    std::cout << "-> Setting efficiency to " << eff << "%" << std::endl;
  }

  float reso_percent = 1.5;
  if( argc>3 ) {
    std::string reso_str(argv[3]);
    reso_percent = atof(reso_str.c_str());
    std::cout << "-> Setting resolution to " << reso_percent << "%" << std::endl;
  }


  //float reso_percent = 2;
  //float reso_percent = 5;

  ZGDrawTools::setStyle();

  std::string suffix = "";
  if( reso_percent!=1.5 ) suffix = std::string(Form("_%.0fpercent", reso_percent));
  std::string dir = cfg.getEventYieldDir() + "_fit_v0_default_shapes_combination_pcorr_lumi_2.3" + suffix;
  std::string fullPath = "../../diphotons/Analysis/macros/" + dir;
  std::string limitsFile( Form( "%s/limits_eff%d.txt", fullPath.c_str(), eff) );

  
  TGraph* gr_obs = new TGraph(0);
  TGraph* gr_exp = new TGraph(0);
  TGraphAsymmErrors* gr_exp_1sigma = new TGraphAsymmErrors(0);
  TGraphAsymmErrors* gr_exp_2sigma = new TGraphAsymmErrors(0);

  ifstream ifs(limitsFile.c_str());
  std::cout << "-> Opened file: " << limitsFile << std::endl;
  int iPoint = 0;
  float lastMass = -1;

  while( ifs.good() ) {

    float m, obs, exp, exp_m1s, exp_m2s, exp_p1s, exp_p2s;
    std::string s_m, s_obs, s_exp, s_exp_m1s, s_exp_m2s, s_exp_p1s, s_exp_p2s;
    ifs >> s_m >> m >> s_obs >> obs >> s_exp >> exp >> s_exp_m1s >> exp_m1s >> s_exp_m2s >> exp_m2s >> s_exp_p1s >>  exp_p1s >> s_exp_p2s >> exp_p2s;
    if( m>=1000. ) continue;
    if( m==lastMass ) continue;
    std::cout << "m: " << m << " obs: " << obs << " exp: " << exp << " exp_m1s: " << exp_m1s << " exp_m2s: " << exp_m2s << " exp_p1s: " << exp_p1s << " exp_p2s: " << exp_p2s << std::endl;

    obs/=2.;
    exp/=2.;
    exp_m1s/=2.;
    exp_m2s/=2.;
    exp_p1s/=2.;
    exp_p2s/=2.;

    gr_obs       ->SetPoint( iPoint, m, obs );
    gr_exp       ->SetPoint( iPoint, m, exp );
    gr_exp_1sigma->SetPoint( iPoint, m, exp );
    gr_exp_2sigma->SetPoint( iPoint, m, exp );

    gr_exp_1sigma->SetPointError( iPoint, 0., 0., exp-exp_m1s, exp_p1s-exp );
    gr_exp_2sigma->SetPointError( iPoint, 0., 0., exp-exp_m2s, exp_p2s-exp );

    iPoint++;
    lastMass = m;

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

  TH2D* h2_axes = new TH2D("axes", "", 10, 350., 950., 10, 0., 9. );
  h2_axes->SetYTitle( "95\% CL UL on #sigma #times BR(A#rightarrowZ#gamma#rightarrowl^{+}l^{-}#gamma) [fb]");
  //h2_axes->SetYTitle( "95\% CL UL [fb]");
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

  c1->SaveAs( Form("%s/limit_eff%d%s.eps", cfg.getEventYieldDir().c_str(), eff, suffix.c_str()) );
  c1->SaveAs( Form("%s/limit_eff%d%s.pdf", cfg.getEventYieldDir().c_str(), eff, suffix.c_str()) );

  return 0;

}
