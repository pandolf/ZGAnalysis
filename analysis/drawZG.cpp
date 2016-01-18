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

  dt.drawPlot( "nVert"     , "nVert"     , "", 35, 0 , 35, "Number of Vertexes" );
  dt.drawPlot( "leptType"     , "leptType"     , "", 5, 9.5 , 14.4, "Lepton PDG ID" );

  dt.drawPlot( "mZg"  , "boss_mass", "", 100, 200., 1200., "M(Z#gamma)", "GeV" );
  dt.drawPlot( "mZg_all"     , "boss_mass", "", 120, 0., 1200., "M(Z#gamma)", "GeV" );

  dt.drawPlot( "mZg_lowptZ"  , "boss_mass", "z_pt<30.", 50, 0., 500., "M(Z#gamma)", "GeV" );
  dt.drawPlot( "mZ_lowptZ"  , "z_mass", "z_pt<30.", 50, 0., 500., "M(Z)", "GeV" );
  dt.drawPlot( "ptGamma_lowptZ"  , "gamma_pt", "z_pt<30.", 50, 0., 500., "Photon p_{T}", "GeV" );
  dt.drawPlot( "ptgOmZg_lowptZ"  , "gamma_pt/boss_mass", "z_pt<30.", 50, 0., 2., "", "" );

  dt.drawPlot( "mZ"   , "z_mass", ""            , 50, 50., 150., "M(l^{+}l^{-})", "GeV" );
  dt.drawPlot( "mZee" , "z_mass", "leptType==11", 50, 50., 150., "M(e^{+}e^{-})", "GeV" );
  dt.drawPlot( "mZmm" , "z_mass", "leptType==13", 50, 50., 150., "M(#mu^{+}#mu^{-})", "GeV" );

  dt.drawPlot( "ptZ"        , "z_pt"    , ""              , 60, 0. , 300., "Z p_{T}"      , "GeV" );
  dt.drawPlot( "ptZ_metCut" , "z_pt"    , "met<50."       , 60, 0. , 300., "Z p_{T}"      , "GeV" );
  dt.drawPlot( "ptZ_bossCut", "z_pt"    , "boss_mass>150. && boss_mass<500.", 60, 0. , 300., "Z p_{T}"      , "GeV" );

  dt.drawPlot( "ptGamma" , "gamma_pt", "", 60, 40., 340., "Photon p_{T}" , "GeV" );
  dt.drawPlot( "etaGamma" , "gamma_eta", "", 50, -3., 3., "Photon #eta" );

  dt.drawPlot( "ptLept0" , "lept0_pt" , "", 60, 25., 325., "Leading Lepton p_{T}" , "GeV" );
  dt.drawPlot( "etaLept0", "lept0_eta", "", 50, -3.,   3., "Leading Lepton #eta" );

  dt.drawPlot( "ptLept1"  , "lept1_pt" , "", 30, 20., 170., "Trailing Lepton p_{T}" , "GeV" );
  dt.drawPlot( "etaLept1" , "lept1_eta", "", 50, -3.,   3., "Trailing Lepton #eta" );

  dt.drawPlot( "met"     , "met"     , "", 60, 0. , 300., "Missing E_{T}", "GeV" );
  dt.drawPlot( "ptgOmZg", "gamma_pt/boss_mass", "", 50, 0., 2., "", "" );

  return 0;

}

