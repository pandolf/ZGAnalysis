#include <iostream>
#include <fstream>

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooBreitWigner.h"

#include "TFile.h"
#include "TH1D.h"
#include "TChain.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLorentzVector.h"

#include "../interface/ZGConfig.h"
#include "../interface/ZGDrawTools.h"
#include "../interface/ZGCommonTools.h"



void fitSingleMass( const std::string& basedir, float mass, const std::string&  width, TGraphErrors* gr_mean, TGraphErrors* gr_width );


int main( int argc, char* argv[] ) {


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


  std::vector<std::string> widths;
  widths.push_back( "0p014" );
  widths.push_back( "1p4" );
  widths.push_back( "5p6" );


  std::string outdir = "genSignalShapes";
  system( Form("mkdir -p %s", outdir.c_str()) );


  TFile* outfile = TFile::Open( Form("%s/genSignalShapes.root", outdir.c_str()), "recreate" );
  outfile->cd();

  for( unsigned iwidth=0; iwidth<widths.size(); ++iwidth ) {

    TGraphErrors* gr_mean  = new TGraphErrors(0);
    TGraphErrors* gr_width = new TGraphErrors(0);

    gr_mean ->SetName( Form("bw_mean_%s" , widths[iwidth].c_str()) );
    gr_width->SetName( Form("bw_width_%s", widths[iwidth].c_str()) );

    for( unsigned imass=0; imass<masses.size(); ++imass ) {

      fitSingleMass( "/pnfs/psi.ch/cms/trivcat/store/user/pandolf/crab/", masses[imass], widths[iwidth], gr_mean, gr_width );

    }

    gr_mean ->Write();
    gr_width->Write();

    delete gr_mean ;
    delete gr_width;

  }


  outfile->Close();

  std::cout << "-> Saved graphs in: " << outfile->GetName() << std::endl;

  return 0;

}





void fitSingleMass( const std::string& basedir, float mass, const std::string&  width, TGraphErrors* gr_mean, TGraphErrors* gr_width ) {


  std::string outdir = "genSignalShapes";
  system( Form("mkdir -p %s", outdir.c_str()) );


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
  int genPart_status[100];
  tree->SetBranchAddress( "genPart_status", genPart_status );



  RooRealVar* x = new RooRealVar("boss_mass", "boss_mass", mass, 0.5*mass, 1.5*mass );

  RooDataSet* data = new RooDataSet( "data", "data", RooArgSet(*x) );



  int nentries = tree->GetEntries();


  for( int iEntry = 0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 25000 == 0 ) std::cout << "  Entry: " << iEntry << " / " << nentries << std::endl;

    tree->GetEntry(iEntry);

    TLorentzVector leptPlus;
    TLorentzVector leptMinus;
    TLorentzVector photon;
    bool foundLeptPlus = false;
    bool foundLeptMinus = false;
    bool foundPhoton = false;
    bool tauEvent = false;

    for( int iPart=0; iPart<ngenPart; ++iPart ) {

      if( genPart_status[iPart]!=1 ) continue;

      if( abs(genPart_pdgId[iPart])==15 ) {
        tauEvent = true;
        break;
      }

      if( genPart_pdgId[iPart]==+11 || genPart_pdgId[iPart]==+13 ) {
        leptMinus.SetPtEtaPhiM( genPart_pt[iPart], genPart_eta[iPart], genPart_phi[iPart], genPart_mass[iPart] );
        foundLeptMinus = true;
      }
      if( genPart_pdgId[iPart]==-11 || genPart_pdgId[iPart]==-13 ) {
        leptPlus.SetPtEtaPhiM( genPart_pt[iPart], genPart_eta[iPart], genPart_phi[iPart], genPart_mass[iPart] );
        foundLeptPlus = true;
      }
      if( genPart_pdgId[iPart]==22 ) {
        photon.SetPtEtaPhiM( genPart_pt[iPart], genPart_eta[iPart], genPart_phi[iPart], genPart_mass[iPart] );
        foundPhoton = true;
      }

    } // for genparts

    if( tauEvent ) continue;
    if( !foundLeptPlus || !foundLeptMinus || !foundPhoton ) continue;


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

    if( photon.Pt()/boss.M()< 40./150. ) continue;

    x->setVal(boss.M());

    data->add(RooArgSet(*x));

  }


  RooRealVar* bw_m = new RooRealVar( "bw_m", "Breit-Wigner Mean" , mass, 0.5*mass, 1.5*mass );
  RooRealVar* bw_w = new RooRealVar( "bw_w", "Breit-Wigner Width", 0.01*mass, 0., 0.5*mass );

  RooBreitWigner* model = new RooBreitWigner( "bw", "Breit-Wigner", *x, *bw_m, *bw_w);

  model->fitTo( *data );

  int npoints = gr_mean->GetN();
  gr_mean ->SetPoint     ( npoints, mass, bw_m->getVal() );
  gr_width->SetPoint     ( npoints, mass, bw_w->getVal() );
  gr_mean ->SetPointError( npoints,   0., bw_m->getError() );
  gr_width->SetPointError( npoints,   0., bw_w->getError() );

  RooPlot* plot = x->frame();
  data->plotOn(plot);
  model->plotOn(plot);

  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();

  plot->Draw();
    
  c1->SaveAs( Form("%s/fit_m%.0f_%s.eps", outdir.c_str(), mass, width.c_str()) );
  c1->SaveAs( Form("%s/fit_m%.0f_%s.pdf", outdir.c_str(), mass, width.c_str()) );

  delete c1;
  delete data;
  delete bw_m;
  delete bw_w;
  delete model;
  delete plot;
  delete x;

}
