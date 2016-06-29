#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>


#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"


#include "interface/ZGSample.h"
#include "interface/ZGConfig.h"
#include "interface/ZGCommonTools.h"


#define zg_cxx
#include "interface/zg.h"


bool doPtWeighting = true;



int main( int argc, char* argv[] ) {


  TChain* tree_clone = new TChain("mt2");
  tree_clone->Add( "$DCAP/pnfs/psi.ch/cms/trivcat/store/user/pandolf/ZGproduction/80X/PostProcessed/mc2016_v2_pu2/skimAndPrune/ZGTo2LG_post_skim.root");
  //tree_clone->Add( "$DCAP/pnfs/psi.ch/cms/trivcat/store/user/pandolf/ZGproduction/80X/PostProcessed/mc2016_v1_pu/skimAndPrune/ZGTo2LG_post_skim.root");
  tree_clone->Add( "$DCAP/pnfs/psi.ch/cms/trivcat/store/user/pandolf/ZGproduction/80X/PostProcessed/mc2016_v2_pu2/skimAndPrune/ZGTo2LG_Pt130_post_skim.root");

  tree_clone->SetBranchStatus("evt_scale1fb", 0);

  ZGTree myTree;
  myTree.Init(tree_clone);


  TChain* tree_read = new TChain("mt2");
  tree_read->Add( "$DCAP/pnfs/psi.ch/cms/trivcat/store/user/pandolf/ZGproduction/80X/PostProcessed/mc2016_v2_pu2/skimAndPrune/ZGTo2LG_post_skim.root");
  //tree_read->Add( "$DCAP/pnfs/psi.ch/cms/trivcat/store/user/pandolf/ZGproduction/80X/PostProcessed/mc2016_v1_pu/skimAndPrune/ZGTo2LG_post_skim.root");
  tree_read->Add( "$DCAP/pnfs/psi.ch/cms/trivcat/store/user/pandolf/ZGproduction/80X/PostProcessed/mc2016_v2_pu2/skimAndPrune/ZGTo2LG_Pt130_post_skim.root");


  TFile* outfile;
  if( doPtWeighting )
    outfile = TFile::Open("ZGTo2LG_post_skim_stitch_ptWeight.root", "recreate");
  else
    outfile = TFile::Open("ZGTo2LG_post_skim_stitch.root", "recreate");
  TTree* newtree = tree_clone->CloneTree(0);

  float evt_scale1fb;
  newtree->Branch("evt_scale1fb", &evt_scale1fb, "evt_scale1fb/F");

  float evt_scale1fb_old;
  tree_read->SetBranchAddress("evt_scale1fb", &evt_scale1fb_old );


  //TF1* f_reweight = new TF1("reweight", "(x>=250.)*(1.1708) + (x<250.)*(exp(0.93627-0.0031*x))", 100., 10000.);
  TF1* f_reweight = new TF1("reweight", "(1.1708)", 100., 10000.);

 
  int nentries = tree_clone->GetEntries();
  
  for( int iEntry=0; iEntry<nentries; ++iEntry ) {
    
    if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;
    
    myTree.GetEntry(iEntry);
    tree_read->GetEntry(iEntry);


    float founPhotonPt = -1.;

    for( unsigned i=0; i<myTree.ngenPart; ++i ) {

      if( myTree.genPart_status[i]!=23 ) continue;
      if( myTree.genPart_pdgId[i]!=23 && myTree.genPart_pdgId[i]!=22 ) continue;
      //if( myTree.genPart_pdgId[i]!=22 ) continue;

      founPhotonPt = myTree.genPart_pt[i];
      break;

    }

    if( founPhotonPt<0. && myTree.evt_id==852 ) continue;
    if( founPhotonPt>0. ) {
      if( founPhotonPt<130. ) {
        if( myTree.evt_id == 852 ) continue;
      } else {
        if( myTree.evt_id == 851 ) continue;
      }
    }


    evt_scale1fb = evt_scale1fb_old;

    if( doPtWeighting && myTree.evt_id==852 ) { // photon pt reweighting
      evt_scale1fb *= 1.17;
      //evt_scale1fb *= f_reweight->Eval( founPhotonPt );
    }

    newtree->Fill();

  }  // for entries

  newtree->Write();
  outfile->Close();

  return 0;

}
