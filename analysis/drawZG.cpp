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


#define DATABLINDING true



int main( int argc, char* argv[] ) {


  if( argc==1 ) {
    std::cout << " USAGE: ./drawZG [cfg] [lumi/shape]" << std::endl;
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
    else if( normType=="onlyMC" || "mcOnly"  ) {
      shapeNorm=false;
      onlyMC=true;
    } else {
      std::cout << "-> Only 'lumi' and 'shape' and 'onlyMC' are supported normTypes." << std::endl;
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

  if( !DATABLINDING || onlyMC ) {
    dt.drawPlot( "mZg_all"     , "boss_mass", "", 120, 0., 1200., "M(Z#gamma)", "GeV" );
    dt.drawPlot( "mZg_lowMass"  , "boss_mass", "", 40, 200., 600., "M(Z#gamma)", "GeV" );
    dt.drawPlot( "mZg"  , "boss_mass", "", 40, 200., 1000., "M(Z#gamma)", "GeV" );
    dt.drawPlot( "mZg_ee_specialBins"  , "boss_mass", "leptType==11", 8, 200., 1400., "M(e^{+}e^{-}#gamma)", "GeV" );
    dt.drawPlot( "mZg_mm_specialBins"  , "boss_mass", "leptType==13", 8, 200., 1400., "M(#mu^{+}#mu^{-}#gamma)", "GeV" );
    dt.drawPlot( "mZg_ee"  , "boss_mass", "leptType==11", 40, 200., 1000., "M(e^{+}e^{-}#gamma)", "GeV" );
    dt.drawPlot( "mZg_mm"  , "boss_mass", "leptType==13", 40, 200., 1000., "M(#mu^{+}#mu^{-}#gamma)", "GeV" );
  }

  dt.drawPlot( "mZ"   , "z_mass", ""            , 50, 50., 150., "M(l^{+}l^{-})", "GeV" );
  dt.drawPlot( "mZee" , "z_mass", "leptType==11", 50, 50., 150., "M(e^{+}e^{-})", "GeV" );
  dt.drawPlot( "mZmm" , "z_mass", "leptType==13", 50, 50., 150., "M(#mu^{+}#mu^{-})", "GeV" );
  dt.drawPlot( "mZee_bossCut" , "z_mass", "leptType==11 && boss_mass>200. && boss_mass<500.", 50, 50., 150., "M(e^{+}e^{-})", "GeV" );
  dt.drawPlot( "mZmm_bossCut" , "z_mass", "leptType==13 && boss_mass>200. && boss_mass<500.", 50, 50., 150., "M(#mu^{+}#mu^{-})", "GeV" );

  dt.drawPlot( "ptZ"        , "z_pt"    , ""              , 30, 0. , 300., "Z p_{T}"      , "GeV" );
  dt.drawPlot( "ptZ_mm"     , "z_pt"    , "leptType==13"  , 30, 0. , 300., "Z p_{T}"      , "GeV" );
  dt.drawPlot( "ptZ_ee"     , "z_pt"    , "leptType==11"  , 30, 0. , 300., "Z p_{T}"      , "GeV" );
  dt.drawPlot( "ptZ_metCut" , "z_pt"    , "met<50."       , 30, 0. , 300., "Z p_{T}"      , "GeV" );
  dt.drawPlot( "ptZ_bossCut", "z_pt"    , "boss_mass>200. && boss_mass<500.", 20, 0. , 300., "Z p_{T}"      , "GeV" );

  dt.drawPlot( "nGamma" , "nGamma", "", 6, 0., 6., "Photon Multiplicity" );
  dt.drawPlot( "ptGamma" , "gamma_pt", "", 60, 40., 340., "Photon p_{T}" , "GeV" );
  dt.drawPlot( "etaGamma" , "gamma_eta", "", 50, -3., 3., "Photon #eta" );
  dt.drawPlot( "isoGamma" , "gamma_iso", "", 50, 0., 10., "Photon Charged Hadron Isolation", "GeV" );

  dt.drawPlot( "ptGamma_mm" , "gamma_pt", "leptType==13", 60, 40., 340., "Photon p_{T}" , "GeV" );
  dt.drawPlot( "ptGamma_ee" , "gamma_pt", "leptType==11", 60, 40., 340., "Photon p_{T}" , "GeV" );

  dt.drawPlot( "ptLept0" , "lept0_pt" , "", 30, 25., 325., "Leading Lepton p_{T}" , "GeV" );
  dt.drawPlot( "etaLept0", "lept0_eta", "", 30, -3.,   3., "Leading Lepton #eta" );

  dt.drawPlot( "ptLept1"  , "lept1_pt" , "", 30, 20., 170., "Trailing Lepton p_{T}" , "GeV" );
  dt.drawPlot( "etaLept1" , "lept1_eta", "", 30, -3.,   3., "Trailing Lepton #eta" );

  dt.drawPlot( "ptMu0" , "lept0_pt" , "leptType==13", 30, 25., 525., "Leading Muon p_{T}" , "GeV" );
  dt.drawPlot( "ptMu1" , "lept1_pt" , "leptType==13", 30, 25., 225., "Trailing Muon p_{T}" , "GeV" );
  dt.drawPlot( "etaMu0" , "lept0_eta", "leptType==13", 30, -3.,   3., "Leading Muon #eta" );
  dt.drawPlot( "etaMu1" , "lept1_eta", "leptType==13", 30, -3.,   3., "Trailing Muon #eta" );

  dt.drawPlot( "ptEle0"  , "lept0_pt" , "leptType==11", 30, 25., 525., "Leading Electron p_{T}" , "GeV" );
  dt.drawPlot( "ptEle1"  , "lept1_pt" , "leptType==11", 30, 25., 225., "Trailing Electron p_{T}" , "GeV" );
  dt.drawPlot( "etaEle0" , "lept0_eta", "leptType==11", 30, -3.,   3., "Leading Electron #eta" );
  dt.drawPlot( "etaEle1" , "lept1_eta", "leptType==11", 30, -3.,   3., "Trailing Electron #eta" );

  //dt.drawPlot( "z_mass_mu370"  , "z_mass"   , "leptType==13 && boss_mass>350. && boss_mass<400.", 25, 50., 130., "M(#mu^{+}#mu^{-})"  , "GeV" );
  //dt.drawPlot( "lept0_pt_mu370", "lept0_pt" , "leptType==13 && boss_mass>350. && boss_mass<400.", 25, 25., 325., "Leading Muon p_{T}" , "GeV" );
  //dt.drawPlot( "lept1_pt_mu370", "lept1_pt" , "leptType==13 && boss_mass>350. && boss_mass<400.", 25, 20., 225., "Trailing Muon p_{T}", "GeV" );
  //dt.drawPlot( "lept0_eta_mu370","lept0_eta", "leptType==13 && boss_mass>350. && boss_mass<400.", 25., -2.5, 2.5,"Leading Muon #eta"  , "" );
  //dt.drawPlot( "lept1_eta_mu370","lept1_eta", "leptType==13 && boss_mass>350. && boss_mass<400.", 25., -2.5, 2.5,"Trailing Muon #eta" , "" );

  dt.drawPlot( "met"     , "met"     , "", 60, 0. , 300., "Missing E_{T}", "GeV" );
  dt.drawPlot( "ptgOmZg", "gamma_pt/boss_mass", "", 25, 0., 1., "p_{T}(#gamma) / M(Z#gamma)", "" );
  dt.drawPlot( "ptgOmZg_bossCut", "gamma_pt/boss_mass", "boss_mass>200. && boss_mass<500.", 25, 0., 1., "p_{T}(#gamma) / M(Z#gamma)", "" );

  return 0;

}

