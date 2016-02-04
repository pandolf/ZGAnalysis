#include <iostream>

#include "RooRealVar.h"
#include "RooCBShape.h"
#include "RooDataSet.h"
#include "RooPlot.h"

#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraphErrors.h"

#include "../interface/ZGDrawTools.h"
#include "../interface/ZGCommonTools.h"



TF1* fitGraph( const std::string& outdir, TGraphErrors* graph, const std::string& axisName );


int main() {


  ZGDrawTools::setStyle();
  
  std::vector<float> masses;
  masses.push_back( 300. );
  masses.push_back( 400. );
  masses.push_back( 500. );
  masses.push_back( 750. );
  masses.push_back( 1000. );
  masses.push_back( 1250. );
  masses.push_back( 1500. );
  masses.push_back( 2000. );
  masses.push_back( 2500. );
  masses.push_back( 3000. );
  masses.push_back( 5000. );


  TFile* outfile = TFile::Open("signalShapeParameters.root", "recreate");
  outfile->cd();

  TGraphErrors* gr_mean = new TGraphErrors(0);
  TGraphErrors* gr_sigma = new TGraphErrors(0);
  TGraphErrors* gr_width = new TGraphErrors(0);
  TGraphErrors* gr_alpha = new TGraphErrors(0);
  TGraphErrors* gr_n = new TGraphErrors(0);

  gr_mean ->SetName("mean");
  gr_sigma->SetName("sigma");
  gr_width->SetName("width");
  gr_alpha->SetName("alpha");
  gr_n    ->SetName("n");


  std::string outdir = "plotsSignalToys";
  system( Form("mkdir -p %s", outdir.c_str() ) );



  for( unsigned i=0; i<masses.size(); ++i ) {

    std::cout << "-> Starting mass: " << masses[i] << std::endl;

    TFile* file  = TFile::Open( Form("toySignal_m%.0f.root", masses[i]) );
    file->cd();

    float thisMass = masses[i];
    RooRealVar x("mass", "mass", thisMass, 0.2*thisMass, 1.5*thisMass );

    // Crystal-Ball
    RooRealVar mean( "mean", "mean", thisMass );
    RooRealVar sigma( "sigma", "sigma", 0.015*thisMass );
    RooRealVar alpha( "alpha", "alpha", 1.2 );
    RooRealVar n( "n", "n", 0.81 );
    RooCBShape cb( "cb", "cb", x, mean, sigma, alpha, n );

    RooDataSet* data = (RooDataSet*)file->Get("data");

    cb.fitTo( *data );

    RooPlot* frame = x.frame();
    data->plotOn(frame);
    cb.plotOn(frame);

    TCanvas* c1 = new TCanvas("c1", "", 600, 600);
    c1->cd();

    frame->Draw();

    c1->SaveAs( Form("%s/fit_m%.0f.eps", outdir.c_str(), thisMass) );
    c1->SaveAs( Form("%s/fit_m%.0f.pdf", outdir.c_str(), thisMass) );

    delete c1;

    gr_mean->SetPoint( i, thisMass, mean.getVal() );
    gr_sigma->SetPoint( i, thisMass, sigma.getVal() );
    gr_width->SetPoint( i, thisMass, sigma.getVal()/mean.getVal() );
    gr_alpha->SetPoint( i, thisMass, alpha.getVal() );
    gr_n->SetPoint( i, thisMass, n.getVal() );

    gr_mean->SetPointError( i, 0., mean.getError() );
    gr_sigma->SetPointError( i, 0., sigma.getError() );
    gr_width->SetPointError( i, 0., sqrt( sigma.getError()*sigma.getError()/(mean.getVal()*mean.getVal()) + sigma.getVal()*sigma.getVal()*mean.getError()*mean.getError()/(mean.getError()*mean.getError()*mean.getError()*mean.getError()) ) );
    gr_alpha->SetPointError( i, 0., alpha.getError() );
    gr_n->SetPointError( i, 0., n.getError() );

  } // for i


  outfile->cd();

  TF1* f1_mean  = fitGraph(outdir, gr_mean , "Gaussian Mean [GeV]");
  TF1* f1_sigma = fitGraph(outdir, gr_sigma, "Gaussian Sigma [GeV]");
  TF1* f1_width = fitGraph(outdir, gr_width, "Gaussian #sigma/#mu");
  TF1* f1_alpha = fitGraph(outdir, gr_alpha, "CB #alpha");
  TF1* f1_n     = fitGraph(outdir, gr_n    , "CB N");

  f1_mean  ->Write();
  f1_sigma ->Write();
  f1_width ->Write();
  f1_alpha ->Write();
  f1_n     ->Write();

  gr_mean  ->Write();
  gr_sigma ->Write();
  gr_alpha ->Write();
  gr_n     ->Write();

  outfile->Close();

  return 0;

}



TF1* fitGraph( const std::string& outdir, TGraphErrors* graph, const std::string& axisName ) {


  float xMin = 200.;
  float xMax = 7000.;

  TF1* f1 = new TF1( Form("f1_%s", graph->GetName()), "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x", xMin, xMax );
  f1->SetLineColor(46);

  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(2);
  graph->SetMarkerColor(46);

  graph->Fit( f1, "QR" );
  TH1D* band = ZGDrawTools::getBand(f1);

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
