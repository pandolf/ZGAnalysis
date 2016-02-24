#include <iostream>

#include "RooRealVar.h"
#include "RooCBShape.h"
#include "RooDataSet.h"
#include "RooPlot.h"

#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TLegend.h"

#include "../interface/ZGConfig.h"
#include "../interface/ZGDrawTools.h"
#include "../interface/ZGCommonTools.h"
//#include "../interface/RooTestCBShape.h"
#include "../interface/RooDoubleCBShape.h"



void fitGraphs( const ZGConfig& cfg, const std::vector<float> masses, const std::string& width, const std::string& outdir, TFile* outfile, const std::string& name, const std::string& sel="" );
TF1* fitGraph( const std::string& outdir, TGraphErrors* graph, const std::string& axisName );
void drawCompare( const ZGConfig& cfg, const std::string& outdir, TFile* file, const std::string& width, const std::string& name, const std::string& axisName, const std::string& cat1, const std::string& name1, const std::string& cat2, const std::string& name2 );


float xMax = 2100.;


int main( int argc, char* argv[] ) {


  if( argc<2 ) {
    std::cout << "USAGE: ./fitSignalShapes [configFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  std::string configFileName(argv[1]);
  ZGConfig cfg(configFileName);


  ZGDrawTools::setStyle();
  
  std::vector<float> masses;
  masses.push_back( 450. );
  masses.push_back( 600. );
  masses.push_back( 750. );
  masses.push_back( 900. );
  masses.push_back( 1500. );
  masses.push_back( 2000. );
  //masses.push_back( 300. );
  //masses.push_back( 400. );
  //masses.push_back( 500. );
  //masses.push_back( 750. );
  //masses.push_back( 1000. );
  //masses.push_back( 1500. );
  //masses.push_back( 2000. );
  //masses.push_back( 2500. );
  //masses.push_back( 3000. );
  //masses.push_back( 5000. );

  std::vector<std::string> widths;
  //widths.push_back( "5p6" );
  widths.push_back( "0p014" );

  for( unsigned iw =0; iw<widths.size(); ++iw ) {

    std::string outdir(Form("%s/signalShapes_w%s", cfg.getEventYieldDir().c_str(), widths[iw].c_str()));
    system( Form("mkdir -p %s", outdir.c_str() ) );

    TFile* outfile = TFile::Open(Form("%s/signalShapeParameters_w%s.root", outdir.c_str(), widths[iw].c_str()), "recreate");
    outfile->cd();

    fitGraphs( cfg, masses, widths[iw], outdir, outfile, "all" );
    fitGraphs( cfg, masses, widths[iw], outdir, outfile, "ee", "leptType==11" );
    fitGraphs( cfg, masses, widths[iw], outdir, outfile, "mm", "leptType==13" );

    drawCompare( cfg, outdir, outfile, widths[iw], "mean"  , "Gaussian Mean [GeV]"   , "ee", "ee#gamma", "mm", "#mu#mu#gamma" );
    drawCompare( cfg, outdir, outfile, widths[iw], "sigma" , "Gaussian #sigma [GeV]" , "ee", "ee#gamma", "mm", "#mu#mu#gamma" );
    drawCompare( cfg, outdir, outfile, widths[iw], "width" , "Gaussian #sigma/#mu"   , "ee", "ee#gamma", "mm", "#mu#mu#gamma" );
    drawCompare( cfg, outdir, outfile, widths[iw], "alpha1", "CB left #alpha"        , "ee", "ee#gamma", "mm", "#mu#mu#gamma" );
    drawCompare( cfg, outdir, outfile, widths[iw], "alpha2", "CB right #alpha"       , "ee", "ee#gamma", "mm", "#mu#mu#gamma" );
    drawCompare( cfg, outdir, outfile, widths[iw], "n1"    , "CB left N"             , "ee", "ee#gamma", "mm", "#mu#mu#gamma" );
    drawCompare( cfg, outdir, outfile, widths[iw], "n2"    , "CB right N"            , "ee", "ee#gamma", "mm", "#mu#mu#gamma" );

    outfile->Close();

  } // for widths

  return 0;

}



void fitGraphs( const ZGConfig& cfg, const std::vector<float> masses, const std::string& width, const std::string& outdir, TFile* outfile, const std::string& name, const std::string& sel ) {


  std::cout << "+++ STARTING: " << name << std::endl;

  TString name_tstr(name);

  std::string plotdir = outdir + "/plots";
  system( Form("mkdir -p %s", plotdir.c_str() ) );

  TGraphErrors* gr_mean   = new TGraphErrors(0);
  TGraphErrors* gr_sigma  = new TGraphErrors(0);
  TGraphErrors* gr_width  = new TGraphErrors(0);
  TGraphErrors* gr_alpha1 = new TGraphErrors(0);
  TGraphErrors* gr_n1     = new TGraphErrors(0);
  TGraphErrors* gr_alpha2 = new TGraphErrors(0);
  TGraphErrors* gr_n2     = new TGraphErrors(0);

  gr_mean  ->SetName(Form("mean_%s" , name.c_str()));
  gr_sigma ->SetName(Form("sigma_%s", name.c_str()));
  gr_width ->SetName(Form("width_%s", name.c_str()));
  gr_alpha1->SetName(Form("alpha1_%s", name.c_str()));
  gr_n1    ->SetName(Form("n1_%s"    , name.c_str()));
  gr_alpha2->SetName(Form("alpha2_%s", name.c_str()));
  gr_n2    ->SetName(Form("n2_%s"    , name.c_str()));




  TFile* file  = TFile::Open( Form("%s/trees.root", cfg.getEventYieldDir().c_str()) );

  outfile->cd();

  for( unsigned i=0; i<masses.size(); ++i ) {

    float thisMass = masses[i];
    std::cout << "-> Starting mass: " << thisMass << std::endl;

    TTree* tree0 = (TTree*)file->Get( Form("XZg_Spin0_ZToLL_W_%s_M_%.0f", width.c_str(), thisMass ));
    if( tree0 == 0 ) continue;
    TTree* tree;
    if( sel=="" ) {
      tree = tree0;
    } else {
      tree = tree0->CopyTree(sel.c_str());
    }

    float maxMassFact = (thisMass>1000. || width=="5p6" ) ? 1.5 : 1.3;
    RooRealVar x("boss_mass", "boss_mass", thisMass, 0.5*thisMass, maxMassFact*thisMass );

    // Crystal-Ball
    RooRealVar mean( "mean", "mean", thisMass, 0.9*thisMass, 1.1*thisMass );
    RooRealVar sigma( "sigma", "sigma", 0.015*thisMass, 0., 0.07*thisMass );
    RooRealVar alpha1( "alpha1", "alpha1", 1.2, 0., 2.5 );
    RooRealVar n1( "n1", "n1", 3., 0., 5. );
    RooRealVar alpha2( "alpha2", "alpha2", 1.2, 0., 2.5 );
    RooRealVar n2( "n2", "n2", 3., 0., 10. );
    RooDoubleCBShape cb( "cb", "cb", x, mean, sigma, alpha1, n1, alpha2, n2 );


    RooDataSet* data = new RooDataSet( "data", "data", RooArgSet(x), RooFit::Import(*tree) );

    cb.fitTo( *data );

    RooPlot* frame = x.frame();
    data->plotOn(frame);
    cb.plotOn(frame);

    TCanvas* c1 = new TCanvas("c1", "", 600, 600);
    c1->cd();

    frame->Draw();

    c1->SaveAs( Form("%s/fit_%s_m%.0f_w%s.eps", plotdir.c_str(), name.c_str(), thisMass, width.c_str()) );
    c1->SaveAs( Form("%s/fit_%s_m%.0f_w%s.pdf", plotdir.c_str(), name.c_str(), thisMass, width.c_str()) );

    c1->SetLogy();

    c1->SaveAs( Form("%s/fit_%s_m%.0f_w%s_log.eps", plotdir.c_str(), name.c_str(), thisMass, width.c_str()) );
    c1->SaveAs( Form("%s/fit_%s_m%.0f_w%s_log.pdf", plotdir.c_str(), name.c_str(), thisMass, width.c_str()) );

    delete c1;

    gr_mean->SetPoint( i, thisMass, mean.getVal() );
    gr_sigma->SetPoint( i, thisMass, sigma.getVal() );
    gr_width->SetPoint( i, thisMass, sigma.getVal()/mean.getVal() );
    gr_alpha1->SetPoint( i, thisMass, alpha1.getVal() );
    gr_n1->SetPoint( i, thisMass, n1.getVal() );
    gr_alpha2->SetPoint( i, thisMass, alpha2.getVal() );
    gr_n2->SetPoint( i, thisMass, n2.getVal() );

    gr_mean->SetPointError( i, 0., mean.getError() );
    gr_sigma->SetPointError( i, 0., sigma.getError() );
    gr_width->SetPointError( i, 0., sigma.getError()/mean.getVal() ); // random
    //gr_width->SetPointError( i, 0., sqrt( sigma.getError()*sigma.getError()/(mean.getVal()*mean.getVal()) + sigma.getVal()*sigma.getVal()*mean.getError()*mean.getError()/(mean.getError()*mean.getError()*mean.getError()*mean.getError()) ) );
    gr_alpha1->SetPointError( i, 0., alpha1.getError() );
    gr_n1->SetPointError( i, 0., n1.getError() );
    gr_alpha2->SetPointError( i, 0., alpha2.getError() );
    gr_n2->SetPointError( i, 0., n2.getError() );

    delete tree;
    delete data;

  } // for i


  TF1* f1_mean   = fitGraph(plotdir, gr_mean  , "Gaussian Mean [GeV]");
  TF1* f1_sigma  = fitGraph(plotdir, gr_sigma , "Gaussian Sigma [GeV]");
  TF1* f1_width  = fitGraph(plotdir, gr_width , "Gaussian #sigma/#mu");
  TF1* f1_alpha1 = fitGraph(plotdir, gr_alpha1, "CB left #alpha");
  TF1* f1_n1     = fitGraph(plotdir, gr_n1    , "CB left N");
  TF1* f1_alpha2 = fitGraph(plotdir, gr_alpha2, "CB right #alpha");
  TF1* f1_n2     = fitGraph(plotdir, gr_n2    , "CB right N");

  f1_mean  ->Write();
  f1_sigma ->Write();
  f1_width ->Write();
  f1_alpha1->Write();
  f1_n1    ->Write();
  f1_alpha2->Write();
  f1_n2    ->Write();

  gr_mean  ->Write();
  gr_sigma ->Write();
  gr_width ->Write();
  gr_alpha1->Write();
  gr_n1    ->Write();
  gr_alpha2->Write();
  gr_n2    ->Write();

}



TF1* fitGraph( const std::string& outdir, TGraphErrors* graph, const std::string& axisName ) {


  //float xMin = 200.;
  //float xMax = 7000.;
  float xMin = 350.;

  TString grName_tstr(graph->GetName());
  std::string formula = "[0] + [1]*x";
  if( grName_tstr.Contains("sigma") && grName_tstr.Contains("mm") && xMax>1100. )
    formula = "[0] + [1]*x + [2]*x*x";
  TF1* f1 = new TF1( Form("f1_%s", graph->GetName()), formula.c_str(), xMin, xMax );
  //TF1* f1 = new TF1( Form("f1_%s", graph->GetName()), "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x", xMin, xMax );
  f1->SetLineColor(46);

  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(2);
  graph->SetMarkerColor(46);

  std::string graphName(graph->GetName());
  TH1D* band;

  //if( grName_tstr.Contains("mm") && ( grName_tstr.Contains("n1") || grName_tstr.Contains("alpha1") ) ) 
  //  f1->SetRange( xMin, 1800. );
  graph->Fit( f1, "QR" );
  f1->SetRange( xMin, xMax );
  band = ZGDrawTools::getBand(f1);

  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();

  float yMax = 0.;
  for( int i=0; i<graph->GetN(); ++i ) {
    Double_t x, y;
    graph->GetPoint(i, x, y);
    if( y>yMax ) yMax = y;
  }

  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., 1.3*yMax );
  h2_axes->SetXTitle( "Mass [GeV]" );
  h2_axes->SetYTitle( axisName.c_str() );
  h2_axes->Draw();

  band->Draw("CE3same");
  graph->Draw("psame");
  
  ZGDrawTools::addLabels(c1, -1., "CMS Simulation");
  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/fit_%s.eps", outdir.c_str(), graph->GetName()) );
  c1->SaveAs( Form("%s/fit_%s.pdf", outdir.c_str(), graph->GetName()) );

  delete c1;
  delete h2_axes;

  return f1;

}



void drawCompare( const ZGConfig& cfg, const std::string& outdir, TFile* file, const std::string& width, const std::string& name, const std::string& axisName, const std::string& cat1, const std::string& name1, const std::string& cat2, const std::string& name2 ){

  std::string plotdir = outdir + "/plots";
  system( Form("mkdir -p %s", plotdir.c_str() ) );

  TGraphErrors* gr1 = (TGraphErrors*)file->Get( Form("%s_%s", name.c_str(), cat1.c_str()) );
  TGraphErrors* gr2 = (TGraphErrors*)file->Get( Form("%s_%s", name.c_str(), cat2.c_str()) );

  float yMax = 0.;
  for( int i=0; i<gr1->GetN(); ++i ) {
    Double_t x, y;
    gr1->GetPoint(i, x, y);
    if( y>yMax ) yMax = y;
    gr2->GetPoint(i, x, y);
    if( y>yMax ) yMax = y;
  }


  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();

  TH2D* h2_axes = new TH2D( "axes", "", 10, 350., xMax, 10, 0., 1.3*yMax );
  h2_axes->SetXTitle( "Generated Mass [GeV]" );
  h2_axes->SetYTitle( axisName.c_str() );
  h2_axes->Draw();

  gr1->SetMarkerStyle(21);
  gr1->SetMarkerSize(1.6);
  gr1->SetMarkerColor(46);
  gr1->SetLineColor(46);
  gr1->GetFunction(Form("f1_%s", gr1->GetName()))->SetLineColor(46);

  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(1.6);
  gr2->SetMarkerColor(38);
  gr2->SetLineColor(38);
  gr2->GetFunction(Form("f1_%s", gr2->GetName()))->SetLineColor(38);

  gr1->Draw("psame");
  gr2->Draw("psame");

  float xMin_leg = 0.6;
  float xMax_leg = 0.9;
  float yMin_leg = 0.2;
  float yMax_leg = 0.38;
  if( name=="width" || name=="sigma" ) {
    xMin_leg = 0.2;
    xMax_leg = 0.5;
    yMin_leg = 0.72;
    yMax_leg = 0.9;
  } 
  if( name=="n2" ) {
    xMin_leg = 0.7;
    xMax_leg = 0.9;
    yMin_leg = 0.7;
    yMax_leg = 0.88;
  } 
  

  std::string widthtext;
  if( width=="0p014" ) widthtext = "W = 0.014%";
  if( width=="1p4"   ) widthtext = "W = 1.4%";
  if( width=="5p6"   ) widthtext = "W = 5.6%";

  TLegend* legend = new TLegend( xMin_leg, yMin_leg, xMax_leg, yMax_leg, widthtext.c_str() );
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  legend->SetFillColor(0);
  legend->AddEntry( gr1, name1.c_str(), "LP" );
  legend->AddEntry( gr2, name2.c_str(), "LP" );
  legend->Draw("same");

  gPad->RedrawAxis();

  ZGDrawTools::addLabels( c1, -1., "CMS Simulation" );

  c1->SaveAs( Form("%s/compare_w%s_%s_%s_vs_%s.eps",  plotdir.c_str(), width.c_str(), name.c_str(), cat1.c_str(), cat2.c_str()) );
  c1->SaveAs( Form("%s/compare_w%s_%s_%s_vs_%s.pdf",  plotdir.c_str(), width.c_str(), name.c_str(), cat1.c_str(), cat2.c_str()) );

  delete c1;
  delete h2_axes;

}
