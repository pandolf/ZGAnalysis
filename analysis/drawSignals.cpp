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


  if( argc==1 ) {
    std::cout << " USAGE: ./drawZG [cfg] [lumi/shape]" << std::endl;
    exit(1);
  }
  

  
  std::string configFileName(argv[1]);
  ZGConfig cfg(configFileName);

  std::string treesFile = cfg.getEventYieldDir() + "/trees.root";

  TFile* file = TFile::Open(treesFile.c_str());

  TTree* tree_M_450 = (TTree*)file->Get("XZg_Spin0_ZToLL_W_0p014_M_450");
  TTree* tree_M_600 = (TTree*)file->Get("XZg_Spin0_ZToLL_W_0p014_M_600");
  TTree* tree_M_750 = (TTree*)file->Get("XZg_Spin0_ZToLL_W_0p014_M_750");
  TTree* tree_M_900 = (TTree*)file->Get("XZg_Spin0_ZToLL_W_0p014_M_900");

  tree_M_450->SetTitle("450 GeV");
  tree_M_600->SetTitle("600 GeV");
  tree_M_750->SetTitle("750 GeV");
  tree_M_900->SetTitle("900 GeV");

  std::string plotsDir = cfg.getEventYieldDir() + "/plotsSignal";
  system( Form("mkdir -p %s", plotsDir.c_str()) );

  std::vector<TTree*> trees_signal;
  trees_signal.push_back( tree_M_450 );
  trees_signal.push_back( tree_M_600 );
  trees_signal.push_back( tree_M_750 );
  trees_signal.push_back( tree_M_900 );

  
  ZGDrawTools dt( plotsDir, cfg.lumi() );
  dt.set_overlay( trees_signal );

  dt.drawPlot( "sig_ptgOmZg", "gamma_pt/boss_mass", "boss_mass>200.", 25, 0., 1., "p_{T}(#gamma) / M(Z#gamma)", "" );

  return 0;

}
