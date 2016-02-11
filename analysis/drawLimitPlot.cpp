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


void drawSingleLimitPlot( const ZGConfig& cfg, const std::string& limitsFile, TF1* f1_eff, float factor, const std::string& cat, const std::string& name, const std::string& axisName );


int main( int argc, char* argv[] ) {


  if( argc<2 ) {
    std::cout << " USAGE: ./drawLimitPlot [cfg] [fitName=\"fit_v0\"]" << std::endl;
    exit(1);
  }
  
  
  std::string configFileName(argv[1]);
  ZGConfig cfg(configFileName);

  std::string fitName = "fit_v0";
  if( argc>2 ) {
    fitName = std::string((argv[2]));
  }


  float reso_percent = 1.5;


  float lumi = 2.3;
  TFile* file_eff = TFile::Open("genAcceptanceTimesEfficiency.root");
  TF1* f1_eff = (TF1*)file_eff->Get("aXe_0p014");
  //float eff = 0.4;
  //float lumi = 20.0;

  //float reso_percent = 2;
  //float reso_percent = 5;

  ZGDrawTools::setStyle();

  std::string suffix = "";
  if( reso_percent!=1.5 ) suffix = std::string(Form("_%.0fpercent", reso_percent));
  //std::string dir( Form("%s_fit_v0_default_shapes_combination_pcorr_lumi_%.1f%s", cfg.getEventYieldDir().c_str(), lumi, suffix.c_str()) );
  //std::string fullPath = "../../diphotons/Analysis/macros/" + dir;
  //std::string limitsFile( Form( "%s/limits_eff%d.txt", fullPath.c_str(), eff) );

  std::string limitsFile( Form( "%s/limits_%s.txt", cfg.getEventYieldDir().c_str(), fitName.c_str() ) );

  drawSingleLimitPlot( cfg, limitsFile, f1_eff, 2.*lumi      , fitName, "Zllgamma", "95\% CL UL on #sigma #times BR(A#rightarrowZ#gamma#rightarrowl^{+}l^{-}#gamma) [fb]" );
  drawSingleLimitPlot( cfg, limitsFile, f1_eff, 2.*lumi*0.033, fitName, "Zgamma"  , "95\% CL UL on #sigma #times BR(A#rightarrowZ#gamma) [fb]" );

  return 0;

}



void drawSingleLimitPlot( const ZGConfig& cfg, const std::string& limitsFile, TF1* f1_eff, float factor, const std::string& cat, const std::string& name, const std::string& axisName ) {

  
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

    float thisEff = f1_eff->Eval(m);
    float conversion = thisEff*factor;
    //float conversion = 2.*eff*lumi;

    obs    /=conversion;
    exp    /=conversion;
    exp_m1s/=conversion;
    exp_m2s/=conversion;
    exp_p1s/=conversion;
    exp_p2s/=conversion;

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

  float yMax = (cfg.lumi()<10.) ? 10. : 2.;
  if( factor<1. ) yMax = 300.;
  if( cat=="fit_ee" || cat=="fit_mm" ) yMax *= 2.;

  TH2D* h2_axes = new TH2D("axes", "", 10, 350., 950., 10, 0., yMax );
  h2_axes->SetYTitle( axisName.c_str() );
  //h2_axes->SetYTitle( "95\% CL UL on #sigma #times BR(A#rightarrowZ#gamma#rightarrowl^{+}l^{-}#gamma) [fb]");
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

  TLegend* legend;
  std::string title = "";
  if( cat=="fit_em" ) {
    title = "ee/#mu#mu Combination";
  } else if( cat=="fit_ee" ) {
    title = "ee Channel";
  } else if( cat=="fit_mm" ) {
    title = "#mu#mu Channel";
  }
  if( title!= "" )
    legend = new TLegend( 0.55, 0.65, 0.9, 0.9, title.c_str() );
  else
    legend = new TLegend( 0.55, 0.7, 0.9, 0.9 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  //legend->AddEntry( gr_exp, "Expected", "L" );
  legend->AddEntry( gr_exp_1sigma, "Expected #pm 1#sigma", "LF" );
  legend->AddEntry( gr_exp_2sigma, "Expected #pm 2#sigma", "LF" );
  //legend->AddEntry( gr_obs, "Observed", "L" );
  legend->Draw("same");


  ZGDrawTools::addLabels( c1, cfg.lumi(), "CMS Simulation");


  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/limit_%s_%s.eps", cfg.getEventYieldDir().c_str(), name.c_str(), cat.c_str()) );
  c1->SaveAs( Form("%s/limit_%s_%s.pdf", cfg.getEventYieldDir().c_str(), name.c_str(), cat.c_str()) );

  delete c1;
  delete h2_axes;

}
