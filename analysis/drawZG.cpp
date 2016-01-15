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

  bool shapeNorm = false;
  if( argc>2 ) {
    std::string normType(argv[2]);
    if( normType=="lumi" ) shapeNorm=false;
    else if( normType=="shape" ) shapeNorm=true;
    else {
      std::cout << "-> Only 'lumi' and 'shape' are supported normTypes." << std::endl;
      exit(17);
    }
  }


  if( shapeNorm )
    std::cout << "-> Using shape normalization." << std::endl;
  else
    std::cout << "-> Using lumi normalization." << std::endl;


  std::string treesFile = cfg.getEventYieldDir() + "/trees.root";

  TFile* file = TFile::Open(treesFile.c_str());

  TTree* tree_data = (TTree*)file->Get("data");
  TTree* tree_zg   = (TTree*)file->Get("zg");
  tree_zg->SetTitle("Z#gamma");

  std::string plotsDir = cfg.getEventYieldDir() + "/plots";
  if( shapeNorm ) plotsDir = plotsDir + "_shape";
  system( Form("mkdir -p %s", plotsDir.c_str()) );

  std::vector<TTree*> trees_mc;
  trees_mc.push_back( tree_zg );


  ZGDrawTools dt(plotsDir, cfg.lumi() );
  dt.set_shapeNorm( shapeNorm );

  dt.set_data( tree_data );
  dt.set_mc( trees_mc );

  dt.set_lumi( cfg.lumi() );

  dt.drawPlot( "mZg"  , "boss_mass", "", 50, 200., 1200., "M(Z#gamma)", "GeV" );

  dt.drawPlot( "mZ"   , "z_mass", ""            , 50, 50., 150., "M(l^{+}l^{-})", "GeV" );
  dt.drawPlot( "mZee" , "z_mass", "leptType==11", 50, 50., 150., "M(e^{+}e^{-})", "GeV" );
  dt.drawPlot( "mZmm" , "z_mass", "leptType==13", 50, 50., 150., "M(#mu^{+}#mu^{-})", "GeV" );

  dt.drawPlot( "ptZ"        , "z_pt"    , ""              , 60, 0. , 300., "Z p_{T}"      , "GeV" );
  dt.drawPlot( "ptZ_metCut" , "z_pt"    , "met<50."       , 60, 0. , 300., "Z p_{T}"      , "GeV" );
  dt.drawPlot( "ptZ_bossCut", "z_pt"    , "boss_mass>200.", 60, 0. , 300., "Z p_{T}"      , "GeV" );

  dt.drawPlot( "ptGamma" , "gamma_pt", "", 60, 40., 340., "Photon p_{T}" , "GeV" );

  dt.drawPlot( "ptLept0" , "lept0_pt", "", 60, 25., 325., "Leading Lepton p_{T}" , "GeV" );
  dt.drawPlot( "ptLept1" , "lept1_pt", "", 30, 20., 170., "Trailing Lepton p_{T}" , "GeV" );

  dt.drawPlot( "met"     , "met"     , "", 60, 0. , 300., "Missing E_{T}", "GeV" );

  return 0;

}

