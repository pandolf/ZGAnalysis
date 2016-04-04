#include <cmath>

#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TEfficiency.h"
#include "TVirtualFitter.h"


#include "../interface/ZGDrawTools.h"


float xMin = 300.;
float xMax = 1000.;

bool use76 = false;


TF1* drawEfficiency( const std::string& outdir, TFile* file, const std::string& suffix1, const std::string& suffix2="", const std::string& legendName1="", const std::string& legendName2="" );
//void drawVsMass( const std::string& outdir, TFile* file, const std::string& suffix="" );
void drawVsMass( const std::string& outdir, TTree* tree, const std::string& suffix="", const std::string& sel="" );
void drawSingleGraph( const std::string& outdir, TGraphErrors* graph, float yMin, float yMax, const std::string& axisName, const std::string& suffix );
void drawResoFromTree( const std::string& outdir, TTree* tree, const std::string& saveName, const std::string& particle, const std::string& var, const std::string& axisName, const std::string& sel="" );



int main() {

  ZGDrawTools::setStyle();

  TFile* file;
  if( use76 )
    file = TFile::Open("genEfficiency76.root");
  else
    file = TFile::Open("genEfficiency.root");

  std::string outdir = "plotsGenEfficiency";
  if( use76 ) outdir += "76";
  system( Form("mkdir -p %s", outdir.c_str()) );

  TFile* outfile;
  if( use76 )
    outfile = TFile::Open( "genEfficiencyFunctions76.root", "recreate" );
  else
    outfile = TFile::Open( "genEfficiencyFunctions.root", "recreate" );
  //TFile* outfile = TFile::Open( Form( "%s/efficiency.root", outdir.c_str() ), "recreate" );

  TF1* f1_all    = drawEfficiency( outdir, file, "all"    );  
  TF1* f1_ee     = drawEfficiency( outdir, file, "ee"     );  
  TF1* f1_mm     = drawEfficiency( outdir, file, "mm"     );  
  TF1* f1_noHLT  = drawEfficiency( outdir, file, "noHLT"  );  
  TF1* f1_noIso  = drawEfficiency( outdir, file, "noIso"  );  

  f1_all   ->Write();
  f1_ee->Write();
  f1_mm->Write();
  f1_noHLT ->Write();
  f1_noIso ->Write();
  outfile->Close();

  drawEfficiency( outdir, file, "ee", "mm"  , "ee#gamma", "#mu#mu#gamma" );
  drawEfficiency( outdir, file, "noIso", "all"  , "Before Photon Isolation", "After Photon Isolation" );  
  drawEfficiency( outdir, file, "noHLT"  , "all", "Before HLT"         , "After HLT" );


  TTree* gentree = (TTree*)file->Get("genTree");

  drawVsMass( outdir, gentree );
  drawVsMass( outdir, gentree, "ee", "leptType==11" );
  drawVsMass( outdir, gentree, "mm", "leptType==13" );


  drawResoFromTree( outdir, gentree, "gammaPt"   , "gamma", "pt"  , "Photon p_{T}" );
  drawResoFromTree( outdir, gentree, "gammaEta"  , "gamma", "eta" , "Photon #eta" );
  drawResoFromTree( outdir, gentree, "gammaPhi"  , "gamma", "phi" , "Photon #phi" );
  drawResoFromTree( outdir, gentree, "gammaMass" , "gamma", "mass", "Photon Mass" );

  drawResoFromTree( outdir, gentree, "lept0Pt"   , "lept0", "pt"  , "Leading Electron p_{T}", "leptType==11" );
  drawResoFromTree( outdir, gentree, "lept0Eta"  , "lept0", "eta" , "Leading Electron #eta" , "leptType==11" );
  drawResoFromTree( outdir, gentree, "lept0Phi"  , "lept0", "phi" , "Leading Electron #phi" , "leptType==11" );
  drawResoFromTree( outdir, gentree, "lept0Mass" , "lept0", "mass", "Leading Electron Mass" , "leptType==11" );

  drawResoFromTree( outdir, gentree, "lept1Pt"   , "lept1", "pt"  , "Trailing Electron p_{T}", "leptType==11" );
  drawResoFromTree( outdir, gentree, "lept1Eta"  , "lept1", "eta" , "Trailing Electron #eta" , "leptType==11" );
  drawResoFromTree( outdir, gentree, "lept1Phi"  , "lept1", "phi" , "Trailing Electron #phi" , "leptType==11" );
  drawResoFromTree( outdir, gentree, "lept1Mass" , "lept1", "mass", "Trailing Electron Mass" , "leptType==11" );

  drawResoFromTree( outdir, gentree, "lept0Pt"   , "lept0", "pt"  , "Leading Muon p_{T}", "leptType==13" );
  drawResoFromTree( outdir, gentree, "lept0Eta"  , "lept0", "eta" , "Leading Muon #eta" , "leptType==13" );
  drawResoFromTree( outdir, gentree, "lept0Phi"  , "lept0", "phi" , "Leading Muon #phi" , "leptType==13" );
  drawResoFromTree( outdir, gentree, "lept0Mass" , "lept0", "mass", "Leading Muon Mass" , "leptType==13" );

  drawResoFromTree( outdir, gentree, "lept1Pt"   , "lept1", "pt"  , "Trailing Muon p_{T}", "leptType==13" );
  drawResoFromTree( outdir, gentree, "lept1Eta"  , "lept1", "eta" , "Trailing Muon #eta" , "leptType==13" );
  drawResoFromTree( outdir, gentree, "lept1Phi"  , "lept1", "phi" , "Trailing Muon #phi" , "leptType==13" );
  drawResoFromTree( outdir, gentree, "lept1Mass" , "lept1", "mass", "Trailing Muon Mass" , "leptType==13" );

  drawResoFromTree( outdir, gentree, "zPt"   , "z", "pt"  , "Z Boson p_{T}" );
  drawResoFromTree( outdir, gentree, "zEta"  , "z", "eta" , "Z Boson #eta" );
  drawResoFromTree( outdir, gentree, "zPhi"  , "z", "phi" , "Z Boson #phi" );
  drawResoFromTree( outdir, gentree, "zMass" , "z", "mass", "Z Boson Mass" );

  drawResoFromTree( outdir, gentree, "bossPt"   , "boss", "pt"  , "Boss p_{T}" );
  drawResoFromTree( outdir, gentree, "bossEta"  , "boss", "eta" , "Boss #eta" );
  drawResoFromTree( outdir, gentree, "bossPhi"  , "boss", "phi" , "Boss #phi" );
  drawResoFromTree( outdir, gentree, "bossMass" , "boss", "mass", "Boss Mass" );

  return 0;

}


TF1* drawEfficiency( const std::string& outdir, TFile* file, const std::string& suffix1, const std::string& suffix2, const std::string& legendName1, const std::string& legendName2 ) {

  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();

  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., 1.001 );
  h2_axes->SetXTitle("Generated Z#gamma Mass [GeV]");
  h2_axes->SetYTitle("Reconstruction Efficiency");
  h2_axes->Draw();

  std::string fullName = "eff";
  if( suffix1 != "" ) fullName = fullName + "_" + suffix1;

  TEfficiency* eff = (TEfficiency*)file->Get(fullName.c_str());
  TGraphAsymmErrors* gr_eff = eff->CreateGraph();
  gr_eff->SetMarkerStyle(20);
  gr_eff->SetMarkerSize(1.3);

  TF1* line = new TF1( Form("line_%s", suffix1.c_str()), "[0]", xMin, xMax );
  //TF1* line = new TF1( Form("line_%s", suffix1.c_str()), "[0] + [1]*x", xMin, xMax );
  line->SetLineColor(46);
  line->SetLineWidth(2);
  gr_eff->Fit( line, "QR0" );

  TH1D* h1_band = ZGDrawTools::getBand( line );
  h1_band->SetFillColor(18);
  if( suffix2=="" ) {
    h1_band->Draw("C E3 same");
    line->Draw("same");
    TPaveText* lineParText = new TPaveText( 0.2, 0.2, 0.55, 0.25, "brNDC" );
    lineParText->SetTextSize(0.038);
    lineParText->SetTextFont(42);
    lineParText->SetFillColor(0);
    lineParText->AddText( Form("Efficiency = %.0f #pm %.0f %%", line->GetParameter(0)*100., line->GetParError(0)*100. ) );
    lineParText->Draw("same");
  }

  gr_eff->Draw("p same");

  if( suffix2!="" ) {

    std::string fullName2 = "eff";
    if( suffix2!="" ) fullName2 = fullName2 + "_" +  suffix2;
    TEfficiency* eff2 = (TEfficiency*)file->Get(fullName2.c_str());
    TGraphAsymmErrors* gr_eff2 = eff2->CreateGraph();
    gr_eff2->SetMarkerStyle(24);
    gr_eff2->SetMarkerSize(1.3);

    gr_eff2->Draw("p same");

    float yMin_legend = ( line->GetParameter(0)>0.5 ) ? 0.2 : 0.75;
    TLegend* legend = new TLegend( 0.2, yMin_legend, 0.55, yMin_legend+0.15 );
    legend->SetFillColor(0);
    legend->SetTextSize(0.038);
    legend->SetTextFont(42);
    legend->AddEntry( gr_eff , legendName1.c_str(), "P" );
    legend->AddEntry( gr_eff2, legendName2.c_str(), "P" );
    legend->Draw("same");

  }

  

  ZGDrawTools::addLabels( c1, -1., "CMS Simulation");

  gPad->RedrawAxis();

  if( suffix2=="" ) {
    c1->SaveAs( Form("%s/gen%s.eps", outdir.c_str(), fullName.c_str()) );
    c1->SaveAs( Form("%s/gen%s.pdf", outdir.c_str(), fullName.c_str()) );
  } else {
    c1->SaveAs( Form("%s/geneff_%s_vs_%s.eps", outdir.c_str(), suffix1.c_str(), suffix2.c_str()) );
    c1->SaveAs( Form("%s/geneff_%s_vs_%s.pdf", outdir.c_str(), suffix1.c_str(), suffix2.c_str()) );
  }

  delete c1;
  delete h2_axes;
  delete h1_band;

  return line;

}



void drawVsMass( const std::string& outdir, TTree* tree, const std::string& suffix, const std::string& sel ) {


  int nBins = 7;
  Double_t bins[nBins+1];
  bins[0] = 300.;
  bins[1] = 350.;
  bins[2] = 400.;
  bins[3] = 450.;
  bins[4] = 500.;
  bins[5] = 600.;
  bins[6] = 700.;
  bins[7] = 1000.;
  //bins[7] = 800.;
  

  std::string suffix2 = suffix;
  if( suffix2 != "" ) suffix2 = "_" + suffix2;

  //std::string getName = "massReso" + suffix2;
  //TH1D* h1_all = (TH1D*)file->Get(getName.c_str());

  TGraphErrors* gr_resp = new TGraphErrors(0);
  TGraphErrors* gr_reso = new TGraphErrors(0);
  gr_resp->SetName( Form("resp%s", suffix2.c_str()) );
  gr_reso->SetName( Form("reso%s", suffix2.c_str()) );


  for( int i=0; i<nBins; i++ ) {

    TCanvas* c1 = new TCanvas( Form("fit_%d", i), "", 600, 600 );
    c1->cd();

    TH1D* h1_reso = new TH1D( Form("reso_%d", i), "", 50, -0.2, 0.2 );
    if( sel!="" )
      tree->Project( h1_reso->GetName(), "(bossReco_mass-bossGen_mass)/bossGen_mass", Form("bossGen_mass>%f && bossGen_mass<%f && %s", bins[i], bins[i+1], sel.c_str()) );
    else
      tree->Project( h1_reso->GetName(), "(bossReco_mass-bossGen_mass)/bossGen_mass", Form("bossGen_mass>%f && bossGen_mass<%f", bins[i], bins[i+1]) );
    //TH1D* h1_reso = (TH1D*)file->Get( Form("reso%s_%d", suffix2.c_str(), i) );
    //if( h1_reso==0 ) break;

    TF1* f1 = new TF1("gaussian", "gaus", -0.2, 0.2);
    f1->SetParameters( h1_reso->Integral(), h1_reso->GetMean(), h1_reso->GetRMS() );
    h1_reso->Fit( f1, "LQR0" );

    int maxiter = 3;
    float nSigma = 1.8;
    for( int iter=0; iter<maxiter; iter++ ) {
      float imean = f1->GetParameter(1);
      float irms = f1->GetParameter(2);
      f1->SetRange( imean-nSigma*irms, imean+nSigma*irms );
      if( iter==maxiter-1 )
        h1_reso->Fit( f1, "LQR" );
      else
        h1_reso->Fit( f1, "LQR0" );
    }

    h1_reso->Draw();
    TPaveText* fitResult = new TPaveText( 0.2, 0.6, 0.5, 0.9, "brNDC" );
    fitResult->SetFillColor(0);
    fitResult->SetTextSize(0.038);
    fitResult->AddText( "Gaussian Fit" );
    fitResult->AddText( Form("#mu = %.2f #pm %.2f %%"   , 100.*f1->GetParameter(1), 100.*f1->GetParError(1)) );
    fitResult->AddText( Form("#sigma = %.2f #pm %.2f %%", 100.*f1->GetParameter(2), 100.*f1->GetParError(2)) );
    fitResult->Draw("same");
    c1->SaveAs( Form("%s/fit%s_%d.eps", outdir.c_str(), suffix2.c_str(), i) );
    c1->SaveAs( Form("%s/fit%s_%d.pdf", outdir.c_str(), suffix2.c_str(), i) );

    float binWidth = bins[i+1]-bins[i];
    float xmass = bins[i] + 0.5*binWidth;
    //float xmass = h1_all->GetBinCenter( i+1 );
    float mean = f1->GetParameter(1);
    float mean_err = f1->GetParError(1);
    gr_resp->SetPoint(i, xmass, mean );
    gr_resp->SetPointError(i, binWidth/2., mean_err );
    //gr_resp->SetPointError(i, binWidth/sqrt(12.), mean_err );

    float sigma = f1->GetParameter(2);
    float sigma_err = f1->GetParError(2);
    gr_reso->SetPoint(i, xmass, sigma );
    gr_reso->SetPointError(i, binWidth/2., sigma_err );
    //gr_reso->SetPointError(i, binWidth/sqrt(12.), sigma_err );
    //gr_reso->SetPointError(i, h1_all->GetBinWidth( i+1 )/sqrt(12.), sqrt( sigma_err*sigma_err/(mean*mean) + sigma*sigma*mean_err*mean_err/(mean*mean*mean*mean) ) );

    delete c1;

  }

  drawSingleGraph( outdir, gr_resp, -0.01, 0.01, "(M_{gen} - M_{reco}) / M_{gen}", suffix );
  drawSingleGraph( outdir, gr_reso,    0., 0.05, "Mass Resolution", suffix );

}


void drawSingleGraph( const std::string& outdir, TGraphErrors* graph, float yMin, float yMax, const std::string& axisName, const std::string& suffix ) {


  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();

  TH2D* h2_axes = new TH2D( "axes", "", 10, xMin, xMax, 10, yMin, yMax );
  h2_axes->SetXTitle( "Generated Z#gamma Mass [GeV]" );
  h2_axes->SetYTitle( axisName.c_str() );
  h2_axes->Draw();

  if( yMin<0. ) {
    TLine* lineHorizontal = new TLine( xMin, 0., xMax, 0. );
    lineHorizontal->Draw("same");
  }
  

  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(1.5);

  TF1* line;
  if( suffix=="ee" ) 
    line = new TF1( "line", "[0]", xMin, xMax );
  else
    line = new TF1( "line", "[0] + [1]*x", xMin, xMax );
  line->SetLineColor(46);
  line->SetLineWidth(2);
  graph->Fit( line, "QR0" );

  TH1D* h1_band = ZGDrawTools::getBand( line );
  h1_band->SetFillColor(18);
  h1_band->Draw("C E3 same");
  line->Draw("same");

  graph->Draw("p same" );

  TPaveText* channel = new TPaveText( 0.2, 0.8, 0.45, 0.9, "brNDC" );
  channel->SetTextSize(0.035);
  channel->SetFillColor(0);
  //channel->SetTextFont(42);
  if( suffix=="ee" ) {
    channel->AddText( "Electron Channel" );
    channel->Draw("same");
  }
  if( suffix=="mm" ) {
    channel->AddText( "Muon Channel" );
    channel->Draw("same");
  }
  
  ZGDrawTools::addLabels( c1, -1., "CMS Simulation" );

  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/gen%s.eps", outdir.c_str(), graph->GetName()) );
  c1->SaveAs( Form("%s/gen%s.pdf", outdir.c_str(), graph->GetName()) );

  delete c1;
  delete h2_axes;
  delete h1_band;

}



void drawResoFromTree( const std::string& outdir, TTree* tree, const std::string& saveName, const std::string& particle, const std::string& var, const std::string& axisName, const std::string& sel ) {


  TH1D* h1_reso = new TH1D("reso", "", 100, -0.2, 0.2 );

  std::string resoVarName(Form("(%sReco_%s-%sGen_%s)/%sGen_%s", particle.c_str(), var.c_str(), particle.c_str(), var.c_str(), particle.c_str(), var.c_str()));
  tree->Project( "reso", resoVarName.c_str(), sel.c_str() );

  TCanvas* c1 = new TCanvas("creso", "", 600, 600);
  c1->cd();

  h1_reso->SetXTitle( Form("%s Resolution", axisName.c_str()) );
  h1_reso->SetYTitle( "Normalized to Unity" );
  h1_reso->DrawNormalized();

  ZGDrawTools::addLabels( c1, -1., "CMS Simulation" );

  c1->SaveAs( Form("%s/reso_%s.eps", outdir.c_str(), saveName.c_str()) );
  c1->SaveAs( Form("%s/reso_%s.pdf", outdir.c_str(), saveName.c_str()) );

  delete h1_reso;
  delete c1;

}
