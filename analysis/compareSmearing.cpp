#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLegend.h"


#include "../interface/ZGDrawTools.h"



class Signal {

 public:
  Signal( float m, const std::string& w) {
    mass = m;
    width = w;
  }

  float mass;
  std::string width;

 private:

};


void drawVar( const std::string& outdir, TTree* tree_yesSmear, TTree* tree_noSmear, const Signal& signal, const std::string& savename, const std::string& name, int nBins, float xMin, float xMax, const std::string& selection="" );

int main() {

  ZGDrawTools::setStyle();

  TFile* file_yesSmear = TFile::Open("EventYields_presel_2016v5/trees.root");
  TFile* file_noSmear  = TFile::Open("EventYields_presel_2016v5_noSmear/trees.root");

  std::vector<Signal> signals;
  signals.push_back(Signal(1000., "0p014"));
  signals.push_back(Signal(2000., "0p014"));
  signals.push_back(Signal( 745., "0p014"));
  signals.push_back(Signal(2000., "5p6"  ));
  signals.push_back(Signal( 500., "5p6"  ));
  //signals.push_back(Signal( 770., "5p6"  ));

  std::string outdir = "compareSmear";
  system( Form("mkdir -p %s", outdir.c_str()) );
  
  for( unsigned i=0; i<signals.size(); ++i ) {

    std::string signalName(Form("GluGluSpin0ToZG_ZToLL_W%s_M%.0f", signals[i].width.c_str(), signals[i].mass));

    TTree* tree_yesSmear = (TTree*)file_yesSmear->Get(signalName.c_str());
    TTree* tree_noSmear  = (TTree*)file_noSmear->Get(signalName.c_str());

    drawVar( outdir, tree_yesSmear, tree_noSmear, signals[i], "gamma_pt","gamma_pt", 100, signals[i].mass*0.2, signals[i].mass*2.);
    drawVar( outdir, tree_yesSmear, tree_noSmear, signals[i], "mZg_ee", "boss_mass", 100, signals[i].mass*0.7, signals[i].mass*1.1, "leptType==11");
    drawVar( outdir, tree_yesSmear, tree_noSmear, signals[i], "mZg_mm", "boss_mass", 100, signals[i].mass*0.7, signals[i].mass*1.1, "leptType==13");
    drawVar( outdir, tree_yesSmear, tree_noSmear, signals[i], "mZee", "z_mass", 80, 50., 130., "leptType==11");
    drawVar( outdir, tree_yesSmear, tree_noSmear, signals[i], "mZmm", "z_mass", 80, 50., 130., "leptType==13");

  }

  return 0;

}



void drawVar( const std::string& outdir, TTree* tree_yesSmear, TTree* tree_noSmear, const Signal& signal, const std::string& savename, const std::string& name, int nBins, float xMin, float xMax, const std::string& selection ) {

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();
 

  float binWidth = (xMax-xMin)/((float)nBins);

  TH1D* h1_yesSmear = new TH1D( "h1_yesSmear", "", nBins, xMin, xMax );
  TH1D* h1_noSmear = new TH1D( "h1_noSmear", "", nBins, xMin, xMax );
  h1_yesSmear->Sumw2();
  h1_noSmear->Sumw2();

  h1_noSmear->SetXTitle( name.c_str() );
  h1_noSmear->SetYTitle( "Events" );

  h1_yesSmear->SetLineWidth(2);
  h1_noSmear->SetLineWidth(2);

  float sf74 = 2.6;
  std::string name74 = name;

  if( selection=="" ) {
    tree_yesSmear->Project("h1_yesSmear", name74.c_str(), Form("weight*%f", sf74));
    tree_noSmear->Project("h1_noSmear", name.c_str(), "weight");
  } else {
    tree_yesSmear->Project("h1_yesSmear", name74.c_str(), Form("weight*%f*(%s)", sf74, selection.c_str()));
    tree_noSmear->Project("h1_noSmear", name.c_str(), Form("weight*(%s)", selection.c_str()));
  }
  h1_noSmear ->SetLineColor(kBlack);
  h1_yesSmear->SetLineColor(kRed);
  h1_noSmear ->DrawNormalized("same");
  h1_yesSmear->DrawNormalized("same");

  ZGDrawTools::addLabels( c1, 2.6 );

  std::string legendTitle( Form( "M=%.0f GeV, %s", signal.mass, signal.width.c_str()) );
  TLegend* legend;
  if( name=="boss_mass" )
    legend = new TLegend( 0.2, 0.72, 0.55, 0.9 );
  else
    legend = new TLegend( 0.65, 0.72, 0.9 , 0.9 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  legend->SetHeader( legendTitle.c_str() );
  legend->AddEntry( h1_noSmear , "Before Smearing", "L" );
  legend->AddEntry( h1_yesSmear, "After Smearing", "L" );
  legend->Draw("same");

  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/%s_M%.0f_W%s.eps"    , outdir.c_str(), savename.c_str(), signal.mass, signal.width.c_str()) );
  c1->SaveAs( Form("%s/%s_M%.0f_W%s.pdf"    , outdir.c_str(), savename.c_str(), signal.mass, signal.width.c_str()) );

  c1->SetLogy();

  c1->SaveAs( Form("%s/%s_M%.0f_W%s_log.eps", outdir.c_str(), savename.c_str(), signal.mass, signal.width.c_str()) );
  c1->SaveAs( Form("%s/%s_M%.0f_W%s_log.pdf", outdir.c_str(), savename.c_str(), signal.mass, signal.width.c_str()) );
  
  delete c1;
  delete h1_yesSmear;
  delete h1_noSmear;

}
