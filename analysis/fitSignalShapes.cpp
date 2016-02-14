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



void fitGraphs( const ZGConfig& cfg, const std::vector<float> masses, const std::string& outdir, TFile* outfile, const std::string& name, const std::string& sel="" );
TF1* fitGraph( const std::string& outdir, TGraphErrors* graph, const std::string& axisName );
void drawCompare( const ZGConfig& cfg, const std::string& outdir, TFile* file, const std::string& name, const std::string& axisName, const std::string& cat1, const std::string& name1, const std::string& cat2, const std::string& name2 );


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
  //masses.push_back( 300. );
  //masses.push_back( 400. );
  //masses.push_back( 500. );
  //masses.push_back( 750. );
  //masses.push_back( 1000. );
  //masses.push_back( 1250. );
  //masses.push_back( 1500. );
  //masses.push_back( 2000. );
  //masses.push_back( 2500. );
  //masses.push_back( 3000. );
  //masses.push_back( 5000. );


  std::string outdir = cfg.getEventYieldDir() + "/signalShapes";
  system( Form("mkdir -p %s", outdir.c_str() ) );

  TFile* outfile = TFile::Open(Form("%s/signalShapeParameters.root", outdir.c_str()), "recreate");
  outfile->cd();

  fitGraphs( cfg, masses, outdir, outfile, "all" );
  fitGraphs( cfg, masses, outdir, outfile, "ee", "leptType==11" );
  fitGraphs( cfg, masses, outdir, outfile, "mm", "leptType==13" );

  drawCompare( cfg, outdir, outfile, "mean" , "Gaussian Mean [GeV]"   , "ee", "ee#gamma", "mm", "#mu#mu#gamma" );
  drawCompare( cfg, outdir, outfile, "sigma", "Gaussian #sigma [GeV]" , "ee", "ee#gamma", "mm", "#mu#mu#gamma" );
  drawCompare( cfg, outdir, outfile, "width", "Gaussian #sigma/#mu"   , "ee", "ee#gamma", "mm", "#mu#mu#gamma" );

  outfile->Close();

  return 0;

}



void fitGraphs( const ZGConfig& cfg, const std::vector<float> masses, const std::string& outdir, TFile* outfile, const std::string& name, const std::string& sel ) {


  std::cout << "+++ STARTING: " << name << std::endl;


  std::string plotdir = outdir + "/plots";
  system( Form("mkdir -p %s", plotdir.c_str() ) );

  TGraphErrors* gr_mean = new TGraphErrors(0);
  TGraphErrors* gr_sigma = new TGraphErrors(0);
  TGraphErrors* gr_width = new TGraphErrors(0);
  TGraphErrors* gr_alpha = new TGraphErrors(0);
  TGraphErrors* gr_n = new TGraphErrors(0);

  gr_mean ->SetName(Form("mean_%s" , name.c_str()));
  gr_sigma->SetName(Form("sigma_%s", name.c_str()));
  gr_width->SetName(Form("width_%s", name.c_str()));
  gr_alpha->SetName(Form("alpha_%s", name.c_str()));
  gr_n    ->SetName(Form("n_%s"    , name.c_str()));




  TFile* file  = TFile::Open( Form("%s/trees.root", cfg.getEventYieldDir().c_str()) );

  outfile->cd();

  for( unsigned i=0; i<masses.size(); ++i ) {

    float thisMass = masses[i];
    std::cout << "-> Starting mass: " << thisMass << std::endl;

    TTree* tree0 = (TTree*)file->Get( Form("XZg_Spin0_ZToLL_W_0p014_M_%.0f", thisMass ));
    TTree* tree;
    if( sel=="" ) {
      tree = tree0;
    } else {
      tree = tree0->CopyTree(sel.c_str());
    }

    RooRealVar x("boss_mass", "boss_mass", thisMass, 0.5*thisMass, 1.2*thisMass );

    // Crystal-Ball
    RooRealVar mean( "mean", "mean", thisMass, 0.9*thisMass, 1.1*thisMass );
    RooRealVar sigma( "sigma", "sigma", 0.015*thisMass, 0., 0.07*thisMass );
    RooRealVar alpha1( "alpha1", "alpha1", 1.2, 0., 2.5 );
    RooRealVar n1( "n1", "n1", 3., 0., 5. );
    RooRealVar alpha2( "alpha2", "alpha2", 1.2, 0., 2.5 );
    RooRealVar n2( "n2", "n2", 3., 0., 5. );
    RooDoubleCBShape cb( "cb", "cb", x, mean, sigma, alpha1, n1, alpha2, n2 );

    //RooRealVar mean( "mean", "mean", thisMass, 0.9*thisMass, 1.1*thisMass );
    //RooRealVar sigma( "sigma", "sigma", 0.015*thisMass, 0., 0.07*thisMass );
    //RooRealVar alpha( "alpha", "alpha", 1.2, 0., 2.5 );
    //RooRealVar n( "n", "n", 3., 0., 5. );
    //RooCBShape cb( "cb", "cb", x, mean, sigma, alpha, n );

    RooDataSet* data = new RooDataSet( "data", "data", RooArgSet(x), RooFit::Import(*tree) );
    //RooDataSet* data = (RooDataSet*)file->Get("data");

    cb.fitTo( *data );

    RooPlot* frame = x.frame();
    data->plotOn(frame);
    cb.plotOn(frame);

    TCanvas* c1 = new TCanvas("c1", "", 600, 600);
    c1->cd();

    frame->Draw();

    c1->SaveAs( Form("%s/fit_%s_m%.0f.eps", plotdir.c_str(), name.c_str(), thisMass) );
    c1->SaveAs( Form("%s/fit_%s_m%.0f.pdf", plotdir.c_str(), name.c_str(), thisMass) );

    c1->SetLogy();

    c1->SaveAs( Form("%s/fit_%s_m%.0f_log.eps", plotdir.c_str(), name.c_str(), thisMass) );
    c1->SaveAs( Form("%s/fit_%s_m%.0f_log.pdf", plotdir.c_str(), name.c_str(), thisMass) );

    delete c1;

    gr_mean->SetPoint( i, thisMass, mean.getVal() );
    gr_sigma->SetPoint( i, thisMass, sigma.getVal() );
    gr_width->SetPoint( i, thisMass, sigma.getVal()/mean.getVal() );
    gr_alpha->SetPoint( i, thisMass, alpha1.getVal() );
    //gr_alpha->SetPoint( i, thisMass, alpha.getVal() );
    gr_n->SetPoint( i, thisMass, n1.getVal() );
    //gr_n->SetPoint( i, thisMass, n.getVal() );

    gr_mean->SetPointError( i, 0., mean.getError() );
    gr_sigma->SetPointError( i, 0., sigma.getError() );
    gr_width->SetPointError( i, 0., sigma.getError()/mean.getVal() ); // random
    //gr_width->SetPointError( i, 0., sqrt( sigma.getError()*sigma.getError()/(mean.getVal()*mean.getVal()) + sigma.getVal()*sigma.getVal()*mean.getError()*mean.getError()/(mean.getError()*mean.getError()*mean.getError()*mean.getError()) ) );
    gr_alpha->SetPointError( i, 0., alpha1.getError() );
    //gr_alpha->SetPointError( i, 0., alpha.getError() );
    gr_n->SetPointError( i, 0., n1.getError() );
    //gr_n->SetPointError( i, 0., n.getError() );

    delete tree;
    delete data;

  } // for i


  TF1* f1_mean  = fitGraph(plotdir, gr_mean , "Gaussian Mean [GeV]");
  TF1* f1_sigma = fitGraph(plotdir, gr_sigma, "Gaussian Sigma [GeV]");
  TF1* f1_width = fitGraph(plotdir, gr_width, "Gaussian #sigma/#mu");
  TF1* f1_alpha = fitGraph(plotdir, gr_alpha, "CB #alpha");
  TF1* f1_n     = fitGraph(plotdir, gr_n    , "CB N");

  f1_mean  ->Write();
  f1_sigma ->Write();
  f1_width ->Write();
  f1_alpha ->Write();
  f1_n     ->Write();

  gr_mean  ->Write();
  gr_sigma ->Write();
  gr_alpha ->Write();
  gr_n     ->Write();
  gr_width ->Write();

}



TF1* fitGraph( const std::string& outdir, TGraphErrors* graph, const std::string& axisName ) {


  //float xMin = 200.;
  //float xMax = 7000.;
  float xMin = 400.;
  float xMax = 1000.;

  TF1* f1 = new TF1( Form("f1_%s", graph->GetName()), "[0] + [1]*x", xMin, xMax );
  //TF1* f1 = new TF1( Form("f1_%s", graph->GetName()), "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x", xMin, xMax );
  f1->SetLineColor(46);

  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(2);
  graph->SetMarkerColor(46);

  std::string graphName(graph->GetName());
  TH1D* band;
  if( graphName!="width" ) {
    graph->Fit( f1, "QR" );
    band = ZGDrawTools::getBand(f1);
  }

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



void drawCompare( const ZGConfig& cfg, const std::string& outdir, TFile* file, const std::string& name, const std::string& axisName, const std::string& cat1, const std::string& name1, const std::string& cat2, const std::string& name2 ){

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

  TH2D* h2_axes = new TH2D( "axes", "", 10, 400., 1000., 10, 0., 1.3*yMax );
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

  TLegend* legend = new TLegend( 0.2, 0.75, 0.55, 0.9 );
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  legend->SetFillColor(0);
  legend->AddEntry( gr1, name1.c_str(), "LP" );
  legend->AddEntry( gr2, name2.c_str(), "LP" );
  legend->Draw("same");

  gPad->RedrawAxis();

  ZGDrawTools::addLabels( c1, -1., "CMS Simulation" );

  c1->SaveAs( Form("%s/compare_%s_%s_vs_%s.eps",  plotdir.c_str(), name.c_str(), cat1.c_str(), cat2.c_str()) );
  c1->SaveAs( Form("%s/compare_%s_%s_vs_%s.pdf",  plotdir.c_str(), name.c_str(), cat1.c_str(), cat2.c_str()) );

  delete c1;
  delete h2_axes;

}
