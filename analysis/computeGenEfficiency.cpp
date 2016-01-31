#include <iostream>
#include <vector>
#include <string>

#include "TGraphErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"




std::vector< TGraphErrors* > getEfficiencyGraphs( const std::string& basedir, std::vector<float> masses );
TGraphErrors* getSingleWidthMasses( const std::string& basedir, std::vector<float> masses, const std::string& width, const std::string& name );



int main() {

  std::vector<float> masses;
  masses.push_back( 300. );
  masses.push_back( 400. );
  masses.push_back( 500. );
  masses.push_back( 650. );
  masses.push_back( 740. );
  masses.push_back( 750. );
  masses.push_back( 760. );
  masses.push_back( 770. );
  masses.push_back( 1000. );
  masses.push_back( 1500. );
  masses.push_back( 2000. );

  std::vector< TGraphErrors* > graphs = getEfficiencyGraphs( "/pnfs/psi.ch/cms/trivcat/store/user/pandolf/crab/", masses );

  return 0;

}



std::vector< TGraphErrors* > getEfficiencyGraphs( const std::string& basedir, std::vector<float> masses ) {

  TGraphErrors* gr_1p4 = getSingleWidthMasses( basedir, masses, "1p4", "W = 1.4\%" );

}

TGraphErrors* getSingleWidthMasses( const std::string& basedir, std::vector<float> masses, const std::string& width, const std::string& name ) {


  TGraphErrors* graph = new TGraphErrors(0);
  graph->SetName( width.c_str() );
  graph->SetTitle( name.c_str() );

  for( unsigned imass = 0; imass<masses.size(); ++imass ) {

    float mass = masses[imass];

    std::string dataset( Form( "GluGluSpin0ToZGamma_ZToLL_W_%d_M_%.0f_TuneCUEP8M1_13TeV_pythia8", width.c_str(), mass ) );
    std::string fileName( Form( "$DCAP/%s/%s/crab_%s/*/0000/genAna_1.root", basedir.c_str(), dataset.c_str(), dataset.c_str() ) );
    std::cout << fileName << std::endl;

  }

  return graph;

}
    
