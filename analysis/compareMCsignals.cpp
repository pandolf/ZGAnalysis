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


void drawVar( const std::string& outdir, TTree* tree_74, TTree* tree_80, const Signal& signal, const std::string& savename, const std::string& name, int nBins, float xMin, float xMax, const std::string& selection="" );

int main() {

  ZGDrawTools::setStyle();

  TFile* file_74 = TFile::Open("~/CMSSW_7_6_3_ZG/src/ZGAnalysis/analysis/EventYields_v0_eth74X/trees.root");
  TFile* file_80 = TFile::Open("EventYields_v0_stitchPt/trees.root");

  std::vector<Signal> signals;
  signals.push_back(Signal(1000., "0p014"));
  signals.push_back(Signal(2000., "0p014"));
  signals.push_back(Signal( 745., "0p014"));
  signals.push_back(Signal(2000., "5p6"  ));
  signals.push_back(Signal( 500., "5p6"  ));
  signals.push_back(Signal( 770., "5p6"  ));

  std::string outdir = "compare74_80";
  system( Form("mkdir -p %s", outdir.c_str()) );
  
  for( unsigned i=0; i<signals.size(); ++i ) {

    std::string signalName(Form("GluGluSpin0ToZG_ZToLL_W%s_M%.0f", signals[i].width.c_str(), signals[i].mass));

    TTree* tree_74;
    if( signals[i].mass==745. && signals[i].width=="0p014" )
      tree_74 = (TTree*)file_74->Get("GluGluSpin0ToZG_ZToLL_W0p014_M750");
    else
      tree_74 = (TTree*)file_74->Get(signalName.c_str());
    TTree* tree_80 = (TTree*)file_80->Get(signalName.c_str());

    drawVar( outdir, tree_74, tree_80, signals[i], "gamma_pt","gamma_pt", 100, signals[i].mass*0.2, signals[i].mass*2.);
    if( signals[i].width=="0p014" ) {
      drawVar( outdir, tree_74, tree_80, signals[i], "mZg_ee", "boss_mass", 100, signals[i].mass*0.7, signals[i].mass*1.2, "leptType==11");
      drawVar( outdir, tree_74, tree_80, signals[i], "mZg_mm", "boss_mass", 100, signals[i].mass*0.7, signals[i].mass*1.2, "leptType==13");
    } else {
      drawVar( outdir, tree_74, tree_80, signals[i], "mZg_ee", "boss_mass", 50, signals[i].mass*0.5, signals[i].mass*1.3, "leptType==11");
      drawVar( outdir, tree_74, tree_80, signals[i], "mZg_mm", "boss_mass", 50, signals[i].mass*0.5, signals[i].mass*1.3, "leptType==13");
    }
    drawVar( outdir, tree_74, tree_80, signals[i], "mZee", "z_mass", 80, 50., 130., "leptType==11");
    drawVar( outdir, tree_74, tree_80, signals[i], "mZmm", "z_mass", 80, 50., 130., "leptType==13");

  }

  return 0;

}



void drawVar( const std::string& outdir, TTree* tree_74, TTree* tree_80, const Signal& signal, const std::string& savename, const std::string& name, int nBins, float xMin, float xMax, const std::string& selection ) {

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();
 

  float binWidth = (xMax-xMin)/((float)nBins);

  TH1D* h1_74 = new TH1D( "h1_74", "", nBins, xMin, xMax );
  TH1D* h1_80 = new TH1D( "h1_80", "", nBins, xMin, xMax );
  h1_74->Sumw2();
  h1_80->Sumw2();

  h1_80->SetXTitle( name.c_str() );
  h1_80->SetYTitle( "Events" );

  h1_74->SetLineWidth(2);
  h1_80->SetLineWidth(2);

  float sf74 = 2.6;
  std::string name74 = name;
  if( signal.mass==745. && signal.width=="0p014" && name=="boss_mass" ) name74 = name + "-5.";

  if( selection=="" ) {
    tree_74->Project("h1_74", name74.c_str(), Form("weight*%f", sf74));
    tree_80->Project("h1_80", name.c_str(), "weight");
  } else {
    tree_74->Project("h1_74", name74.c_str(), Form("weight*%f*(%s)", sf74, selection.c_str()));
    tree_80->Project("h1_80", name.c_str(), Form("weight*(%s)", selection.c_str()));
  }
  h1_80->SetLineColor(kRed);
  h1_80->DrawNormalized("same");
  h1_74->DrawNormalized("same");

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
  legend->AddEntry( h1_74, "74X MC", "L" );
  legend->AddEntry( h1_80, "80X MC", "L" );
  legend->Draw("same");

  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/%s_M%.0f_W%s_74vs80.eps"    , outdir.c_str(), savename.c_str(), signal.mass, signal.width.c_str()) );
  c1->SaveAs( Form("%s/%s_M%.0f_W%s_74vs80.pdf"    , outdir.c_str(), savename.c_str(), signal.mass, signal.width.c_str()) );

  c1->SetLogy();

  c1->SaveAs( Form("%s/%s_M%.0f_W%s_74vs80_log.eps", outdir.c_str(), savename.c_str(), signal.mass, signal.width.c_str()) );
  c1->SaveAs( Form("%s/%s_M%.0f_W%s_74vs80_log.pdf", outdir.c_str(), savename.c_str(), signal.mass, signal.width.c_str()) );
  
  delete c1;
  delete h1_74;
  delete h1_80;

}
