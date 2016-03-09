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


void drawSingleLimitPlot( const ZGConfig& cfg, const std::string& limitsFile, TF1* f1_eff, const std::string& cat, const std::string& width, const std::string& name, bool onlyExpected=false );


int main( int argc, char* argv[] ) {


  if( argc<2 ) {
    std::cout << " USAGE: ./drawLimitPlot [cfg] [fitName=\"fit_v0\"] [width=\"0p014\"]" << std::endl;
    exit(1);
  }
  
  
  std::string configFileName(argv[1]);
  ZGConfig cfg(configFileName);

  std::string fitName = "fit_v0";
  if( argc>2 ) {
    fitName = std::string((argv[2]));
  }

  std::string width = "0p014";
  if( argc>3 ) {
    width = std::string((argv[3]));
  }




  //TFile* file_eff = TFile::Open("genAcceptanceTimesEfficiency.root");
  //TF1* f1_eff = (TF1*)file_eff->Get("f1_gr_0p014_times_line_all");

  TFile* file_eff = TFile::Open(Form("%s/signalEfficiency_w%s.root", cfg.getEventYieldDir().c_str(), width.c_str()));
  TF1* f1_eff = (TF1*)file_eff->Get("f1_eff_all");

  ZGDrawTools::setStyle();


  std::string limitsFile( Form( "%s/limits_w%s_%s.txt", cfg.getEventYieldDir().c_str(), width.c_str(), fitName.c_str() ) );

  drawSingleLimitPlot( cfg, limitsFile, f1_eff, fitName, width, "Zllgamma", true );
  drawSingleLimitPlot( cfg, limitsFile, f1_eff, fitName, width, "Zllgamma", false );
  drawSingleLimitPlot( cfg, limitsFile, f1_eff, fitName, width, "Zgamma"  , true   );
  drawSingleLimitPlot( cfg, limitsFile, f1_eff, fitName, width, "Zgamma"  , false   );

  return 0;

}





void drawSingleLimitPlot( const ZGConfig& cfg, const std::string& limitsFile, TF1* f1_eff, const std::string& cat, const std::string& width, const std::string& name, bool onlyExpected ) {


  std::string axisName;
  float factor;
  if( name=="Zgamma" ) {
    axisName = "95\% CL UL on #sigma #times BR(A#rightarrowZ#gamma) [fb]";
    factor = 2.*cfg.lumi()*0.0337;
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
  int iPointExp = 0;
  int iPointObs = 0;
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


    if( obs>0. ) {
      gr_obs       ->SetPoint( iPointObs, m, obs );
      iPointObs++;
    }

    bool okForExp = ((m==348.) || (m==2000.) || (m==748.)
                  || (m>=400. && m<900.  && (int(m) % 50 == 0) && m!=750.) 
                  //|| (m>=700. && m<1100. && (int(m) % 100 == 0))
                  || (m>=900. && (int(m-100) % 200 == 0))  );
    if( okForExp ) {
      gr_exp       ->SetPoint( iPointExp, m, exp );
      gr_exp_1sigma->SetPoint( iPointExp, m, exp );
      gr_exp_2sigma->SetPoint( iPointExp, m, exp );
      gr_exp_1sigma->SetPointError( iPointExp, 0., 0., exp-exp_m1s, exp_p1s-exp );
      gr_exp_2sigma->SetPointError( iPointExp, 0., 0., exp-exp_m2s, exp_p2s-exp );
      iPointExp++;
    }

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
  if( width=="5p6" ) yMax *= 2.;

  TH2D* h2_axes = new TH2D("axes", "", 10, 350., 1000., 10, 0., yMax );
  h2_axes->SetYTitle( axisName.c_str() );
  //h2_axes->SetYTitle( "95\% CL UL on #sigma #times BR(A#rightarrowZ#gamma#rightarrowl^{+}l^{-}#gamma) [fb]");
  h2_axes->SetXTitle( "Resonance Mass [GeV]");
  h2_axes->Draw();


  gr_exp_2sigma->Draw("E3 same");
  gr_exp_1sigma->Draw("E3 same");
  gr_exp       ->Draw("L  same");
  if( !onlyExpected )
    gr_obs       ->Draw("L  same");


  gr_exp_1sigma->SetLineWidth(2);
  gr_exp_1sigma->SetLineStyle(2);
  gr_exp_2sigma->SetLineWidth(2);
  gr_exp_2sigma->SetLineStyle(2);

  TLegend* legend;
  std::string title = "";
  if( width=="0p014" ) title = "W = 0.014%";
  if( width=="5p6"   ) title = "W = 5.6%";
  if( onlyExpected )
    legend = new TLegend( 0.55, 0.65, 0.9, 0.9 );
  else
    legend = new TLegend( 0.55, 0.6 , 0.9, 0.9 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  if( title!= "" )
    legend->SetHeader(title.c_str());
  //legend->AddEntry( gr_exp, "Expected", "L" );
  if( !onlyExpected )
    legend->AddEntry( gr_obs, "Observed", "L" );
  legend->AddEntry( gr_exp_1sigma, "Expected #pm 1#sigma", "LF" );
  legend->AddEntry( gr_exp_2sigma, "Expected #pm 2#sigma", "LF" );
  legend->Draw("same");


  if( onlyExpected )
    ZGDrawTools::addLabels( c1, cfg.lumi(), "CMS Simulation");
  else
    ZGDrawTools::addLabels( c1, cfg.lumi(), "CMS Preliminary");


  gPad->RedrawAxis();

  if( onlyExpected ) {
    c1->SaveAs( Form("%s/limitExp_%s_w%s_%s.eps", cfg.getEventYieldDir().c_str(), name.c_str(), width.c_str(), cat.c_str()) );
    c1->SaveAs( Form("%s/limitExp_%s_w%s_%s.pdf", cfg.getEventYieldDir().c_str(), name.c_str(), width.c_str(), cat.c_str()) );
  } else {
    c1->SaveAs( Form("%s/limit_%s_w%s_%s.eps", cfg.getEventYieldDir().c_str(), name.c_str(), width.c_str(), cat.c_str()) );
    c1->SaveAs( Form("%s/limit_%s_w%s_%s.pdf", cfg.getEventYieldDir().c_str(), name.c_str(), width.c_str(), cat.c_str()) );
  }

  c1->Clear();

  TH2D* h2_axes_long = new TH2D("axes_long", "", 10, 350., 2000., 10, 0., yMax );
  h2_axes_long->SetYTitle( axisName.c_str() );
  h2_axes_long->SetXTitle( "Resonance Mass [GeV]");
  h2_axes_long->Draw();


  gr_exp_2sigma->Draw("E3 same");
  gr_exp_1sigma->Draw("E3 same");
  gr_exp       ->Draw("L  same");
  if( !onlyExpected )
    gr_obs       ->Draw("L  same");


  legend->Draw("same");


  if( onlyExpected )
    ZGDrawTools::addLabels( c1, cfg.lumi(), "CMS Simulation");
  else
    ZGDrawTools::addLabels( c1, cfg.lumi(), "CMS Preliminary");

  gPad->RedrawAxis();

  if( onlyExpected ) {
    c1->SaveAs( Form("%s/limitExp_%s_w%s_%s_long.eps", cfg.getEventYieldDir().c_str(), name.c_str(), width.c_str(), cat.c_str()) );
    c1->SaveAs( Form("%s/limitExp_%s_w%s_%s_long.pdf", cfg.getEventYieldDir().c_str(), name.c_str(), width.c_str(), cat.c_str()) );
  } else {
    c1->SaveAs( Form("%s/limit_%s_w%s_%s_long.eps", cfg.getEventYieldDir().c_str(), name.c_str(), width.c_str(), cat.c_str()) );
    c1->SaveAs( Form("%s/limit_%s_w%s_%s_long.pdf", cfg.getEventYieldDir().c_str(), name.c_str(), width.c_str(), cat.c_str()) );
  }


  delete c1;
  delete h2_axes;
  delete h2_axes_long;

}
