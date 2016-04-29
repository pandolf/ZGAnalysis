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
    std::cout << " USAGE: ./drawQ [cfg] [lumi/shape]" << std::endl;
    exit(1);
  }
  

  
  std::string configFileName(argv[1]);
  ZGConfig cfg(configFileName);

  bool shapeNorm = false;
  bool onlyMC = false;
  if( argc>2 ) {
    std::string normType(argv[2]);
    if( normType=="lumi" ) shapeNorm=false;
    else if( normType=="shape" ) shapeNorm=true;
    else if( normType=="onlyMC" ) {
      shapeNorm=false;
      onlyMC=true;
    } else {
      std::cout << "-> Only 'lumi' and 'shape' are supported normTypes." << std::endl;
      exit(17);
    }
  }


  if( shapeNorm )
    std::cout << "-> Using shape normalization." << std::endl;
  else
    std::cout << "-> Using lumi normalization." << std::endl;

  if( onlyMC )
    std::cout << "-> Plotting only MC." << std::endl;


  std::string treesFile = cfg.getEventYieldDir() + "/trees.root";

  TFile* file = TFile::Open(treesFile.c_str());
  if( file==0 ) {
    std::cout << "-> File " << treesFile << " doesn't seem to exist!" << std::endl;
    exit(1);
  }

  TTree* tree_data = (TTree*)file->Get("data");
  TTree* tree_zg   = (TTree*)file->Get("zg");
  tree_zg->SetTitle("Z#gamma");
  TTree* tree_dy   = (TTree*)file->Get("dy");
  tree_dy->SetTitle("Z+jets");
  TTree* tree_top   = (TTree*)file->Get("top");
  if( tree_top )
    tree_top->SetTitle("Top");

  std::string plotsDir = cfg.getEventYieldDir() + "/plots";
  if( shapeNorm ) plotsDir = plotsDir + "_shape";
  if( onlyMC ) plotsDir = plotsDir + "_MConly";
  system( Form("mkdir -p %s", plotsDir.c_str()) );

  std::vector<TTree*> trees_mc;
  trees_mc.push_back( tree_zg );
  trees_mc.push_back( tree_dy );
  //trees_mc.push_back( tree_top );


  ZGDrawTools dt(plotsDir, cfg.lumi() );
  dt.set_shapeNorm( shapeNorm );

  if( !onlyMC ) {
    dt.set_data( tree_data );
  }
  dt.set_mc( trees_mc );

  dt.set_lumi( cfg.lumi() );

  dt.drawPlot( "nVert"     , "nVert"     , "", 35, 0 , 35, "Number of Vertexes" );
  dt.drawPlot( "nVert_noPU"     , "nVert"     , "1./puWeight", 35, 0 , 35, "Number of Vertexes" );
  dt.drawPlot( "leptType"     , "leptType"     , "", 5, 9.5 , 14.5, "Lepton PDG ID" );
  dt.drawPlot( "leptType_bossCut"     , "leptType"     , "boss_mass>200. && boss_mass<500.", 5, 9.5 , 14.5, "Lepton PDG ID" );

  dt.drawPlot( "mQgamma_all" , "qgamma_mass", "", 60, 0., 1200., "M(q#gamma)", "GeV" );
  dt.drawPlot( "mQZ_all"     , "qZ_mass"    , "", 60, 0., 1200., "M(qZ)", "GeV" );
  dt.drawPlot( "mQgamma" , "qgamma_mass", "", 40, 150., 950., "M(q#gamma)", "GeV" );
  dt.drawPlot( "mQZ_ee"  , "qZ_mass", "leptType==11", 40, 200., 1000., "M(e^{+}e^{-}q)", "GeV" );
  dt.drawPlot( "mQZ_mm"  , "qZ_mass", "leptType==13", 40, 200., 1000., "M(#mu^{+}#mu^{-}q)", "GeV" );


  return 0;

}

