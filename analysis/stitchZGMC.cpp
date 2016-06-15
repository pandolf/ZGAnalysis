#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>


#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TRandom3.h"


#include "interface/ZGSample.h"
#include "interface/ZGConfig.h"
#include "interface/ZGCommonTools.h"


#define zg_cxx
#include "interface/zg.h"




int main( int argc, char* argv[] ) {


  TChain* tree = new TChain("mt2");
  tree->Add( "$DCAP/pnfs/psi.ch/cms/trivcat/store/user/pandolf/ZGproduction/80X/PostProcessed/mc2016_v2/skimAndPrune/ZGTo2LG_post_skim.root");
  tree->Add( "$DCAP/pnfs/psi.ch/cms/trivcat/store/user/pandolf/ZGproduction/80X/PostProcessed/mc2016_v2/skimAndPrune/ZGTo2LG_Pt130_post_skim.root");


  ZGTree myTree;
  myTree.Init(tree);


  TFile* outfile = TFile::Open("ZGTo2LG_post_skim_stitch.root", "recreate");
  TTree* newtree = tree->CloneTree(0);


 
  int nentries = tree->GetEntries();
  
  for( int iEntry=0; iEntry<nentries; ++iEntry ) {
    
    if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;
    
    myTree.GetEntry(iEntry);


    float founPhotonPt = -1.;

    for( unsigned i=0; i<myTree.ngenPart; ++i ) {

      if( myTree.genPart_status[i]!=23 ) continue;
      if( myTree.genPart_pdgId[i]!=22 ) continue;

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

    newtree->Fill();

  }  // for entries

  newtree->Write();
  outfile->Close();

  return 0;

}
