#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TH1D.h"
#include "TH2D.h"

#include "../interface/ZGDrawTools.h"
#include "../interface/ZGConfig.h"





int main( int argc, char* argv[] ) {


  //if( argc==1 ) {
  //  std::cout << " USAGE: ./drawSignals [cfg]" << std::endl;
  //  exit(1);
  //}
  
  //std::string configFileName(argv[1]);
  //ZGConfig cfg(configFileName);


  TFile* file_M_300  = TFile::Open("genAcceptance/0p014/tree_m300.root");
  TFile* file_M_500  = TFile::Open("genAcceptance/0p014/tree_m500.root");
  TFile* file_M_650  = TFile::Open("genAcceptance/0p014/tree_m650.root");
  TFile* file_M_750  = TFile::Open("genAcceptance/0p014/tree_m750.root");
  TFile* file_M_1000 = TFile::Open("genAcceptance/0p014/tree_m1000.root");
  TFile* file_M_1500 = TFile::Open("genAcceptance/0p014/tree_m1500.root");
  TFile* file_M_2000 = TFile::Open("genAcceptance/0p014/tree_m2000.root");

  TTree* tree_M_300  = (TTree*)file_M_300 ->Get("genTree_M_300");
  TTree* tree_M_500  = (TTree*)file_M_500 ->Get("genTree_M_500");
  TTree* tree_M_650  = (TTree*)file_M_650 ->Get("genTree_M_650");
  TTree* tree_M_750  = (TTree*)file_M_750 ->Get("genTree_M_750");
  TTree* tree_M_1000 = (TTree*)file_M_1000->Get("genTree_M_1000");
  TTree* tree_M_1500 = (TTree*)file_M_1500->Get("genTree_M_1500");
  TTree* tree_M_2000 = (TTree*)file_M_2000->Get("genTree_M_2000");

  tree_M_300 ->SetTitle("300 GeV");
  tree_M_500 ->SetTitle("500 GeV");
  tree_M_650 ->SetTitle("650 GeV");
  tree_M_750 ->SetTitle("750 GeV");
  tree_M_1000->SetTitle("1000 GeV");
  tree_M_1500->SetTitle("1500 GeV");
  tree_M_2000->SetTitle("2000 GeV");

  std::string plotsDir = "plotsGenSignal";
  //std::string plotsDir = cfg.getEventYieldDir() + "/plotsGenSignal";
  system( Form("mkdir -p %s", plotsDir.c_str()) );

  std::vector<TTree*> trees_signal;
  trees_signal.push_back( tree_M_300 );
  trees_signal.push_back( tree_M_500 );
  //trees_signal.push_back( tree_M_650 );
  trees_signal.push_back( tree_M_750 );
  trees_signal.push_back( tree_M_1000 );
  trees_signal.push_back( tree_M_1500 );
  trees_signal.push_back( tree_M_2000 );

  
  ZGDrawTools dt( plotsDir, 2.6 );
  dt.set_overlay( trees_signal );
  dt.set_yAxisScale( 1.1 );
  dt.set_addOverflow( false );

  dt.drawPlot( "gensig_ptgOmZg", "gamma_pt/boss_mass", "boss_mass>200.", 25, 0., 1., "p_{T}(#gamma) / M(Z#gamma)", "" );
  dt.drawPlot( "gensig_deltaRLept", "deltaR_lept", "boss_mass>200.", 50, 0., 2., "Generated #DeltaR (leptons)", "" );

  return 0;

}
