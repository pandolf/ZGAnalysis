#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "TGraphErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TLegend.h"

#include "../interface/ZGDrawTools.h"




std::vector< TGraphErrors* > getEfficiencyGraphs( const std::string& basedir, std::vector<float> masses );
TGraphErrors* getSingleWidthMasses( const std::string& basedir, std::vector<float> masses, const std::string& width, const std::string& name );



float xMin = 200.;
float xMax = 7100.;


int main() {

  std::vector<float> masses;
  masses.push_back( 300. );
  masses.push_back( 400. );
  masses.push_back( 500. );
  masses.push_back( 650. );
  masses.push_back( 740. );
  masses.push_back( 745. );
  masses.push_back( 750. );
  masses.push_back( 755. );
  masses.push_back( 760. );
  masses.push_back( 765. );
  masses.push_back( 770. );
  masses.push_back( 1000. );
  masses.push_back( 1500. );
  masses.push_back( 2000. );
  masses.push_back( 3000. );
  masses.push_back( 4000. );
  masses.push_back( 5000. );
  masses.push_back( 6000. );
  masses.push_back( 7000. );

  std::vector< TGraphErrors* > graphs = getEfficiencyGraphs( "/pnfs/psi.ch/cms/trivcat/store/user/pandolf/crab/", masses );

  ZGDrawTools::setStyle();

  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();
  
  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., 1.0001 );
  h2_axes->SetXTitle( "Mass [GeV]" );
  h2_axes->SetYTitle( "Generator Acceptance" );
  h2_axes->Draw();

  TLegend* legend = new TLegend( 0.2, 0.7, 0.5, 0.91 );
  legend->SetFillColor( 0 );
  legend->SetTextFont( 42 );
  legend->SetTextSize( 0.038 );

  std::vector<int> colors;
  colors.push_back(46);
  colors.push_back(42);
  colors.push_back(38);

  for( unsigned i=0; i<graphs.size(); ++i ) {

    TF1* f1 = new TF1( Form("f1_%d", i), "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x", xMin, xMax );
    f1->SetLineColor( colors[i] );
    graphs[i]->Fit( f1, "QR" );

    graphs[i]->SetMarkerSize( 1.3 );
    graphs[i]->SetMarkerStyle( 20+i );
    graphs[i]->SetMarkerColor( colors[i] );
    graphs[i]->Draw(" p same " );

    legend->AddEntry( graphs[i], graphs[i]->GetTitle(), "P" );

  }

  legend->Draw("same");
  
  ZGDrawTools::addLabels( c1, -1., "CMS Simulation");

  gPad->RedrawAxis(); 

  c1->SaveAs( "genAcceptance.eps" );
  c1->SaveAs( "genAcceptance.pdf" );
  
  c1->SetLogx();
  h2_axes->GetXaxis()->SetNoExponent();
  h2_axes->GetXaxis()->SetMoreLogLabels();
  c1->SaveAs( "genAcceptance_logx.eps" );
  c1->SaveAs( "genAcceptance_logx.pdf" );
  

//TFile* file = TFile::Open("prova.root", "recreate" );
//file->cd();
//for( unsigned i=0; i<graphs.size(); ++i ) 
//  graphs[i]->Write();
//file->Close();

  return 0;

}




std::vector< TGraphErrors* > getEfficiencyGraphs( const std::string& basedir, std::vector<float> masses ) {

  std::vector<TGraphErrors*> graphs;
  graphs.push_back( getSingleWidthMasses( basedir, masses, "0p014", "W = 0.014\%" ) );
  graphs.push_back( getSingleWidthMasses( basedir, masses, "1p4", "W = 1.4\%" ) );
  graphs.push_back( getSingleWidthMasses( basedir, masses, "5p6", "W = 5.6\%" ) );

  return graphs;

}




TGraphErrors* getSingleWidthMasses( const std::string& basedir, std::vector<float> masses, const std::string& width, const std::string& name ) {


  TGraphErrors* graph = new TGraphErrors(0);
  graph->SetName( Form("gr_%s", width.c_str()) );
  graph->SetTitle( name.c_str() );

  for( unsigned imass = 0; imass<masses.size(); ++imass ) {

    float mass = masses[imass];

    if( width=="1p4" ) {
      if( mass==1000. || mass==2000. ) {
        std::cout << " Skipping..." << std::endl;
        continue;
      }
    }

    std::string dataset( Form( "GluGluSpin0ToZGamma_ZToLL_W_%s_M_%.0f_TuneCUEP8M1_13TeV_pythia8", width.c_str(), mass ) );
    std::cout << "-> Starting: " << dataset << std::endl;

    system( Form("ls %s/%s/crab_%s/*/0000/genAna_1.root >> toBeAdded.txt", basedir.c_str(), dataset.c_str(), dataset.c_str() ) );
    TChain* tree = new TChain("mt2");
    ifstream ifs("toBeAdded.txt");
    while( ifs.good() ) {
      std::string fileName;
      ifs >> fileName;
      TString fileName_tstr(fileName);
      if( !fileName_tstr.Contains("pnfs") ) continue;
      tree->Add(Form("$DCAP/%s/mt2", fileName.c_str()) );
    }
    system( "rm toBeAdded.txt" );



    int ngenPart;
    tree->SetBranchAddress( "ngenPart", &ngenPart );
    float genPart_pt[100];
    tree->SetBranchAddress( "genPart_pt", genPart_pt );
    float genPart_eta[100];
    tree->SetBranchAddress( "genPart_eta", genPart_eta );
    float genPart_phi[100];
    tree->SetBranchAddress( "genPart_phi", genPart_phi );
    float genPart_mass[100];
    tree->SetBranchAddress( "genPart_mass", genPart_mass );
    int genPart_pdgId[100];
    tree->SetBranchAddress( "genPart_pdgId", genPart_pdgId );


    int nentries = tree->GetEntries();

    float eff_num = 0.;
    float eff_denom = 0.;


    for( int iEntry = 0; iEntry<nentries; ++iEntry ) {

      if( iEntry % 25000 == 0 ) std::cout << "  Entry: " << iEntry << " / " << nentries << std::endl;

      tree->GetEntry(iEntry);

      eff_denom += 1.;

      TLorentzVector leptPlus;
      TLorentzVector leptMinus;
      TLorentzVector photon;

      for( int iPart=0; iPart<ngenPart; ++iPart ) {

        if( genPart_pdgId[iPart]==+11 || genPart_pdgId[iPart]==+13 )
          leptMinus.SetPtEtaPhiM( genPart_pt[iPart], genPart_eta[iPart], genPart_phi[iPart], genPart_mass[iPart] );
        if( genPart_pdgId[iPart]==-11 || genPart_pdgId[iPart]==-13 )
          leptPlus.SetPtEtaPhiM( genPart_pt[iPart], genPart_eta[iPart], genPart_phi[iPart], genPart_mass[iPart] );
        if( genPart_pdgId[iPart]==22 )
          photon.SetPtEtaPhiM( genPart_pt[iPart], genPart_eta[iPart], genPart_phi[iPart], genPart_mass[iPart] );

      } // for genparts

      if( photon.Pt() < 40. ) continue;
      float ptMax = TMath::Max( leptPlus.Pt(), leptMinus.Pt() );
      float ptMin = TMath::Min( leptPlus.Pt(), leptMinus.Pt() );
      if( ptMax<25. ) continue;
      if( ptMin<20. ) continue;
      if( fabs( photon.Eta() ) > 2.5 ) continue;
      if( fabs( photon.Eta())>1.44 && fabs( photon.Eta())<1.57 ) continue;
      if( fabs( leptPlus.Eta() ) > 2.4 ) continue;
      if( fabs( leptMinus.Eta() ) > 2.4 ) continue;
      if( photon.DeltaR( leptPlus  ) < 0.4 ) continue;
      if( photon.DeltaR( leptMinus ) < 0.4 ) continue;

      TLorentzVector zBoson = leptPlus + leptMinus;
      if( zBoson.M() < 50. ) continue;

      TLorentzVector boss = zBoson + photon;
      if( boss.M() < 200. ) continue;

      eff_num += 1.;

    } // for entries

    float eff = (eff_denom>0. ) ? eff_num/eff_denom : -1.;
    if( eff<0. ) continue;
    graph->SetPoint( imass, mass, eff );
    graph->SetPointError( imass, 0., sqrt( eff*(1.-eff)/eff_denom ) );

    delete tree;


  } // for masses


  return graph;

}
    
