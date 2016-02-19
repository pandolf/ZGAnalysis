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


void drawSingleLimitPlot( const ZGConfig& cfg, const std::string& limitsFile, TF1* f1_eff, const std::string& cat, const std::string& name );


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


  //TFile* file_eff = TFile::Open("genAcceptanceTimesEfficiency.root");
  //TF1* f1_eff = (TF1*)file_eff->Get("f1_gr_0p014_times_line_all");

  TFile* file_eff = TFile::Open(Form("%s/signalEfficiency.root", cfg.getEventYieldDir().c_str()));
  TF1* f1_eff = (TF1*)file_eff->Get("f1_0p014_all");

  ZGDrawTools::setStyle();

  std::string suffix = "";
  if( reso_percent!=1.5 ) suffix = std::string(Form("_%.0fpercent", reso_percent));
  //std::string dir( Form("%s_fit_v0_default_shapes_combination_pcorr_lumi_%.1f%s", cfg.getEventYieldDir().c_str(), lumi, suffix.c_str()) );
  //std::string fullPath = "../../diphotons/Analysis/macros/" + dir;
  //std::string limitsFile( Form( "%s/limits_eff%d.txt", fullPath.c_str(), eff) );

  std::string limitsFile( Form( "%s/limits_%s.txt", cfg.getEventYieldDir().c_str(), fitName.c_str() ) );

  drawSingleLimitPlot( cfg, limitsFile, f1_eff, fitName, "Zllgamma" );
  drawSingleLimitPlot( cfg, limitsFile, f1_eff, fitName, "Zgamma"   );

  return 0;

}





void drawSingleLimitPlot( const ZGConfig& cfg, const std::string& limitsFile, TF1* f1_eff, const std::string& cat, const std::string& name ) {


  std::string axisName;
  float factor;
  if( name=="Zgamma" ) {
    axisName = "95\% CL UL on #sigma #times BR(A#rightarrowZ#gamma) [fb]";
    factor = 2.*cfg.lumi()*0.033;
  } else if( name=="Zllgamma" ) {
    std::string leptType = "l";
    if( cat=="fit_em" || cat=="fit_v0" ) leptType="l";
    if( cat=="fit_ee" ) leptType="e";
    if( cat=="fit_mm" ) leptType="#mu";
    axisName = std::string(Form("95\%% CL UL on #sigma #times BR(A#rightarrowZ#gamma#rightarrow %s^{+}%s^{-}#gamma) [fb]", leptType.c_str(), leptType.c_str()));
    factor = 2.*cfg.lumi();
  } else {
    std::cout << "UNKNOWN NAME! (" << name << ")" << std::endl;
    return;
  }


  
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
    TString m_tstr(s_m);
    if( m_tstr.BeginsWith("#") ) continue;
    if( m==lastMass ) continue;
    if( m==350. ) continue;

    float thisEff = f1_eff->Eval(m);
    float conversion = thisEff*factor;
    //float conversion = 2.*eff*lumi;

    obs    /=conversion;
    exp    /=conversion;
    exp_m1s/=conversion;
    exp_m2s/=conversion;
    exp_p1s/=conversion;
    exp_p2s/=conversion;
    std::cout << "m: " << m << " obs: " << obs << " exp: " << exp << " exp_m1s: " << exp_m1s << " exp_m2s: " << exp_m2s << " exp_p1s: " << exp_p1s << " exp_p2s: " << exp_p2s << std::endl;

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

  TH2D* h2_axes = new TH2D("axes", "", 10, 350., 1000., 10, 0., yMax );
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
  title = "Narrow Resonance";
  //if( cat=="fit_em" ) {
  //  //title = "ee/#mu#mu Combination";
  //} else if( cat=="fit_ee" ) {
  //  title = "Narrow Resonance (e^{+}e^{-}#gamma)";
  //} else if( cat=="fit_mm" ) {
  //  title = "Narrow Resonance (#mu^{+}#mu^{-}#gamma)";
  //  //title = "#mu#mu Channel";
  //}
  legend = new TLegend( 0.55, 0.65, 0.9, 0.9 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  if( title!= "" )
    legend->SetHeader(title.c_str());
  //legend->AddEntry( gr_exp, "Expected", "L" );
  legend->AddEntry( gr_exp_1sigma, "Expected #pm 1#sigma", "LF" );
  legend->AddEntry( gr_exp_2sigma, "Expected #pm 2#sigma", "LF" );
  //legend->AddEntry( gr_obs, "Observed", "L" );
  legend->Draw("same");


  ZGDrawTools::addLabels( c1, cfg.lumi(), "CMS Simulation");


  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/limit_%s_%s.eps", cfg.getEventYieldDir().c_str(), name.c_str(), cat.c_str()) );
  c1->SaveAs( Form("%s/limit_%s_%s.pdf", cfg.getEventYieldDir().c_str(), name.c_str(), cat.c_str()) );

  c1->Clear();

  TH2D* h2_axes_long = new TH2D("axes_long", "", 10, 350., 2000., 10, 0., yMax );
  h2_axes_long->SetYTitle( axisName.c_str() );
  h2_axes_long->SetXTitle( "Resonance Mass [GeV]");
  h2_axes_long->Draw();


  gr_exp_2sigma->Draw("E3 same");
  gr_exp_1sigma->Draw("E3 same");
  gr_exp       ->Draw("L  same");
  //gr_obs       ->Draw("L  same");


  legend->Draw("same");

  ZGDrawTools::addLabels( c1, cfg.lumi(), "CMS Simulation");

  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/limit_%s_%s_long.eps", cfg.getEventYieldDir().c_str(), name.c_str(), cat.c_str()) );
  c1->SaveAs( Form("%s/limit_%s_%s_long.pdf", cfg.getEventYieldDir().c_str(), name.c_str(), cat.c_str()) );

  delete c1;
  delete h2_axes;
  delete h2_axes_long;

}
