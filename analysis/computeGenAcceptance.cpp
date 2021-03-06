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




std::vector< TGraphErrors* > getEfficiencyGraphs( const std::string& outdir, const std::string& basedir, std::vector<float> masses, bool applyRunningCut=true );
TGraphErrors* getSingleWidthMasses( const std::string& outdir, const std::string& basedir, std::vector<float> masses, const std::string& width, const std::string& name, bool applyCut=true );
void drawRelativeEfficiency( std::vector<TGraphErrors*> graphs_denom, std::vector<TGraphErrors*> graphs_num );



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


  std::string outputdir = "genAcceptance";
  system( Form("mkdir -p %s", outputdir.c_str()) );

  std::vector< TGraphErrors* > graphs_noCut = getEfficiencyGraphs( outputdir, "/pnfs/psi.ch/cms/trivcat/store/user/pandolf/crab/", masses, false );
  std::vector< TGraphErrors* > graphs       = getEfficiencyGraphs( outputdir, "/pnfs/psi.ch/cms/trivcat/store/user/pandolf/crab/", masses, true  );

  ZGDrawTools::setStyle();

  drawRelativeEfficiency( graphs_noCut, graphs );

  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();
  
  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, 6500, 10, 0., 1.0001 );
  h2_axes->SetXTitle( "Mass [GeV]" );
  h2_axes->SetYTitle( "Generator Acceptance" );
  h2_axes->Draw();

  TLegend* legend = new TLegend( 0.2, 0.2, 0.5, 0.41 );
  legend->SetFillColor( 0 );
  legend->SetTextFont( 42 );
  legend->SetTextSize( 0.038 );

  std::vector<int> colors;
  colors.push_back(46);
  colors.push_back(42);
  colors.push_back(38);


  std::vector<TF1*> functions;

  for( unsigned i=0; i<graphs.size(); ++i ) {

    TF1* f1 = new TF1( Form("f1_%s", graphs[i]->GetName()), "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x", xMin, xMax );
    f1->SetLineColor( colors[i] );
    graphs[i]->Fit( f1, "QR" );

    graphs[i]->SetMarkerSize( 1.3 );
    graphs[i]->SetMarkerStyle( 20+i );
    graphs[i]->SetMarkerColor( colors[i] );
    graphs[i]->Draw(" p same " );

    legend->AddEntry( graphs[i], graphs[i]->GetTitle(), "P" );

    functions.push_back(f1);

  }

  legend->Draw("same");
  
  ZGDrawTools::addLabels( c1, -1., "CMS Simulation");

  gPad->RedrawAxis(); 

  c1->SaveAs( Form("%s/genAcceptance.eps", outputdir.c_str()) );
  c1->SaveAs( Form("%s/genAcceptance.pdf", outputdir.c_str()) );
  
  c1->SetLogx();
  h2_axes->GetXaxis()->SetNoExponent();
  h2_axes->GetXaxis()->SetMoreLogLabels();
  c1->SaveAs( Form("%s/genAcceptance_logx.eps", outputdir.c_str()) );
  c1->SaveAs( Form("%s/genAcceptance_logx.pdf", outputdir.c_str()) );

  TFile* outfile = TFile::Open( Form("%s/genAcceptance.root", outputdir.c_str()), "recreate");
  outfile->cd();
  for( unsigned i=0; i<functions.size(); ++i )
    functions[i]->Write();
  outfile->Close();
  

  std::cout << "-> Saved stuff in : " << outfile->GetName() << std::endl;

//TFile* file = TFile::Open("prova.root", "recreate" );
//file->cd();
//for( unsigned i=0; i<graphs.size(); ++i ) 
//  graphs[i]->Write();
//file->Close();

  return 0;

}




std::vector< TGraphErrors* > getEfficiencyGraphs( const std::string& outdir, const std::string& basedir, std::vector<float> masses, bool applyRunningCut ) {

  std::vector<TGraphErrors*> graphs;
  graphs.push_back( getSingleWidthMasses( outdir, basedir, masses, "0p014", "W = 0.014\%", applyRunningCut ) );
  graphs.push_back( getSingleWidthMasses( outdir, basedir, masses, "1p4"  , "W = 1.4\%"  , applyRunningCut ) );
  graphs.push_back( getSingleWidthMasses( outdir, basedir, masses, "5p6"  , "W = 5.6\%"  , applyRunningCut ) );

  return graphs;

}




TGraphErrors* getSingleWidthMasses( const std::string& outdir, const std::string& basedir, std::vector<float> masses, const std::string& width, const std::string& name, bool applyCut ) {


  std::string outputdir(Form("%s/%s", outdir.c_str(), width.c_str()));
  system( Form("mkdir -p %s", outputdir.c_str()) );

  TGraphErrors* graph = new TGraphErrors(0);
  if( applyCut )
    graph->SetName( Form("gr_%s", width.c_str()) );
  else
    graph->SetName( Form("gr_%s_noCut", width.c_str()) );
  graph->SetTitle( name.c_str() );

  for( unsigned imass = 0; imass<masses.size(); ++imass ) {

    float mass = masses[imass];

    //if( width=="1p4" ) {
    //  if( mass==1000. || mass==2000. ) {
    //    std::cout << " Skipping..." << std::endl;
    //    continue;
    //  }
    //}

    TFile* outfile = TFile::Open( Form("%s/tree_m%.0f.root", outputdir.c_str(), mass), "recreate" );
    TTree* outTree = new TTree(Form("genTree_M_%.0f", mass), "");

    float weight;
    outTree->Branch( "weight", &weight, "weight/F");

    int leptType;
    outTree->Branch( "leptType", &leptType, "leptType/I");

    float lept0_pt;
    outTree->Branch( "lept0_pt", &lept0_pt, "lept0_pt/F" );
    float lept0_eta;
    outTree->Branch( "lept0_eta", &lept0_eta, "lept0_eta/F" );
    float lept0_phi;
    outTree->Branch( "lept0_phi", &lept0_phi, "lept0_phi/F" );

    float lept1_pt;
    outTree->Branch( "lept1_pt", &lept1_pt, "lept1_pt/F" );
    float lept1_eta;
    outTree->Branch( "lept1_eta", &lept1_eta, "lept1_eta/F" );
    float lept1_phi;
    outTree->Branch( "lept1_phi", &lept1_phi, "lept1_phi/F" );

    float deltaR_lept;
    outTree->Branch( "deltaR_lept", &deltaR_lept, "deltaR_lept/F" );

    float gamma_pt;
    outTree->Branch( "gamma_pt", &gamma_pt, "gamma_pt/F" );
    float gamma_eta;
    outTree->Branch( "gamma_eta", &gamma_eta, "gamma_eta/F" );
    float gamma_phi;
    outTree->Branch( "gamma_phi", &gamma_phi, "gamma_phi/F" );
    float gamma_iso;
    outTree->Branch( "gamma_iso", &gamma_iso, "gamma_iso/F" );

    float z_pt;
    outTree->Branch( "z_pt", &z_pt, "z_pt/F" );
    float z_eta;
    outTree->Branch( "z_eta", &z_eta, "z_eta/F" );
    float z_phi;
    outTree->Branch( "z_phi", &z_phi, "z_phi/F" );
    float z_mass;
    outTree->Branch( "z_mass", &z_mass, "z_mass/F" );

    float boss_pt;
    outTree->Branch( "boss_pt", &boss_pt, "boss_pt/F" );
    float boss_eta;
    outTree->Branch( "boss_eta", &boss_eta, "boss_eta/F" );
    float boss_phi;
    outTree->Branch( "boss_phi", &boss_phi, "boss_phi/F" );
    float boss_mass;
    outTree->Branch( "boss_mass", &boss_mass, "boss_mass/F" );



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


    int nentries = tree->GetEntries();

    float eff_denom = 0.;
    float eff_num = 0.;


    for( int iEntry = 0; iEntry<nentries; ++iEntry ) {

      if( iEntry % 25000 == 0 ) std::cout << "  Entry: " << iEntry << " / " << nentries << std::endl;

      tree->GetEntry(iEntry);

      weight = 1.;

      TLorentzVector leptPlus;
      TLorentzVector leptMinus;
      TLorentzVector photon;
      bool foundLeptPlus = false;
      bool foundLeptMinus = false;
      bool foundPhoton = false;
      bool tauEvent = false;
      leptType = 0;

      for( int iPart=0; iPart<ngenPart; ++iPart ) {

        if( genPart_status[iPart]!=1 ) continue;

        if( abs(genPart_pdgId[iPart])==15 ) {
          tauEvent = true;
          break;
        }

        if( genPart_pdgId[iPart]==+11 || genPart_pdgId[iPart]==+13 ) {
          leptMinus.SetPtEtaPhiM( genPart_pt[iPart], genPart_eta[iPart], genPart_phi[iPart], genPart_mass[iPart] );
          foundLeptMinus = true;
          leptType = abs(genPart_pdgId[iPart]);
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

      eff_denom += 1.;


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
      if( zBoson.M() < 50. || zBoson.M()>130. ) continue;

      TLorentzVector boss = zBoson + photon;
      if( boss.M() < 200. ) continue;

      if( applyCut ) {
        if( photon.Pt()/boss.M()< 40./150. ) continue;
        //if( fabs(photon.Eta())>1.44 ) continue;
      }

      if( leptPlus.Pt() > leptMinus.Pt() ) {

        lept0_pt  = leptPlus.Pt();
        lept0_eta = leptPlus.Eta();
        lept0_phi = leptPlus.Phi();

        lept1_pt  = leptMinus.Pt();
        lept1_eta = leptMinus.Eta();
        lept1_phi = leptMinus.Phi();

      } else {

        lept0_pt  = leptMinus.Pt();
        lept0_eta = leptMinus.Eta();
        lept0_phi = leptMinus.Phi();

        lept1_pt  = leptPlus.Pt();
        lept1_eta = leptPlus.Eta();
        lept1_phi = leptPlus.Phi();

      }

      deltaR_lept = leptPlus.DeltaR(leptMinus);

      gamma_pt  = photon.Pt();
      gamma_eta = photon.Eta();
      gamma_phi = photon.Phi();

      z_pt   = zBoson.Pt();
      z_eta  = zBoson.Eta();
      z_phi  = zBoson.Phi();
      z_mass = zBoson.M();

      boss_pt   = boss.Pt();
      boss_eta  = boss.Eta();
      boss_phi  = boss.Phi();
      boss_mass = boss.M();

      outTree->Fill();

      eff_num += 1.;


    } // for entries

    outfile->cd();
    outTree->Write();
    outfile->Close();
    std::cout << "-> Wrote tree to file: " << outfile->GetName() << std::endl;

    float eff = (eff_denom>0. ) ? eff_num/eff_denom : -1.;
    if( eff<0. ) continue;

    graph->SetPoint( imass, mass, eff );
    graph->SetPointError( imass, 0., sqrt( eff*(1.-eff)/eff_denom ) );

    delete tree;


  } // for masses


  return graph;

}
    


void drawRelativeEfficiency( std::vector<TGraphErrors*> graphs_denom, std::vector<TGraphErrors*> graphs_num ) {

  TCanvas* c1 = new TCanvas( "c2", "", 600, 600 );
  c1->cd();

  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., 1.0001 );
  h2_axes->SetXTitle( "Mass [GeV]" );
  h2_axes->SetYTitle( "Selection Efficiency" );
  h2_axes->Draw();

  TLegend* legend = new TLegend( 0.2, 0.2, 0.5, 0.41 );
  legend->SetFillColor( 0 );
  legend->SetTextFont( 42 );
  legend->SetTextSize( 0.038 );

  std::vector<int> colors;
  colors.push_back(46);
  colors.push_back(42);
  colors.push_back(38);



  for( unsigned i=0; i<graphs_denom.size(); ++i ) {

    TGraphErrors* gr_ratio = new TGraphErrors(0);
    gr_ratio->SetName( Form("%s_cutEff", graphs_num[i]->GetName()) );
    gr_ratio->SetTitle( graphs_num[i]->GetTitle() );

    for( int iPoint=0; iPoint<graphs_denom[i]->GetN(); ++iPoint ) {

      Double_t x,effnum;
      graphs_num[i]->GetPoint(iPoint, x, effnum);
      Double_t effnum_err = graphs_num[i]->GetErrorY(iPoint);

      Double_t effdenom;
      graphs_denom[i]->GetPoint(iPoint, x, effdenom);
      Double_t effdenom_err = graphs_denom[i]->GetErrorY(iPoint);

      Double_t eff = effnum/effdenom;
      Double_t eff_err = sqrt( effnum_err*effnum_err/(effdenom*effdenom) + effnum*effnum*effdenom_err*effdenom_err/(effdenom*effdenom*effdenom*effdenom) );
      gr_ratio->SetPoint( iPoint, x, eff );
      gr_ratio->SetPointError( iPoint, 0., eff_err );

    }

    TF1* f1 = new TF1( Form("f2_%s", gr_ratio->GetName()), "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x", xMin, xMax );
    f1->SetLineColor( colors[i] );
    gr_ratio->Fit( f1, "QR" );

    gr_ratio->SetMarkerSize( 1.3 );
    gr_ratio->SetMarkerStyle( 20+i );
    gr_ratio->SetMarkerColor( colors[i] );
    gr_ratio->Draw(" p same " );

    legend->AddEntry( gr_ratio, gr_ratio->GetTitle(), "P" );

  }

  legend->Draw("same");
  
  ZGDrawTools::addLabels( c1, -1., "CMS Simulation");

  gPad->RedrawAxis(); 

  std::string outputdir = "genAcceptance";
  system( Form("mkdir -p %s", outputdir.c_str()) );

  c1->SaveAs( Form("%s/cutEffSignal.eps", outputdir.c_str()) );
  c1->SaveAs( Form("%s/cutEffSignal.pdf", outputdir.c_str()) );
  
  c1->SetLogx();
  h2_axes->GetXaxis()->SetNoExponent();
  h2_axes->GetXaxis()->SetMoreLogLabels();
  c1->SaveAs( Form("%s/cutEffSignal_logx.eps", outputdir.c_str()) );
  c1->SaveAs( Form("%s/cutEffSignal_logx.pdf", outputdir.c_str()) );

  delete c1;
  delete h2_axes;

}

