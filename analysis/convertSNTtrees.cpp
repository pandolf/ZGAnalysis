#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TLorentzVector.h"





bool loose_skim = false;


void convertFiles( const std::string& samples, const std::string& expr, int id, bool skim=false, bool prune=false );



int main( int argc, char* argv[] ) {


  if( argc<2 ) {
    std::cout << "USAGE: ./convertSNTtrees [samples]" << std::endl;
    exit(11);
  }

  std::string samples(argv[1]);
  //std::string samples = "Run2015_25ns74X";

  bool doSkim = true;
  bool doPrune = false;

  convertFiles(samples, "DoubleEG"      , 4, doSkim, doPrune);
  convertFiles(samples, "DoubleMuon"    , 5, doSkim, doPrune);
  convertFiles(samples, "SinglePhoton"  , 7, doSkim, doPrune);
  convertFiles(samples, "SingleMuon"    , 8, doSkim, doPrune);
  convertFiles(samples, "SingleElectron", 9, doSkim, doPrune);

  //convertFiles(samples, "ZGTo2LG",      851, doSkim, doPrune);
  convertFiles(samples, "GluGluSpin0ToZG_ZToLL_W0p014_M1000", 14251, doSkim, doPrune );
  convertFiles(samples, "GluGluSpin0ToZG_ZToLL_W0p014_M1250", 14250, doSkim, doPrune );
  convertFiles(samples, "GluGluSpin0ToZG_ZToLL_W0p014_M1500", 14249, doSkim, doPrune );
  convertFiles(samples, "GluGluSpin0ToZG_ZToLL_W0p014_M1750", 14248, doSkim, doPrune );
  convertFiles(samples, "GluGluSpin0ToZG_ZToLL_W0p014_M2000", 14247, doSkim, doPrune );
  convertFiles(samples, "GluGluSpin0ToZG_ZToLL_W0p014_M500" , 14236, doSkim, doPrune );
  convertFiles(samples, "GluGluSpin0ToZG_ZToLL_W0p014_M750" , 14229, doSkim, doPrune );
  convertFiles(samples, "GluGluSpin0ToZG_ZToLL_W5p6_M2000"  , 14214, doSkim, doPrune );
  convertFiles(samples, "GluGluSpin0ToZG_ZToLL_W5p6_M500"   , 14210, doSkim, doPrune );
  convertFiles(samples, "GluGluSpin0ToZG_ZToLL_W5p6_M750"   , 14205, doSkim, doPrune );
  convertFiles(samples, "GluGluSpin0ToZG_ZToLL_W5p6_M770"   , 14201, doSkim, doPrune );

  convertFiles(samples, "XZg_Spin0_ZToLL_W_0p014_M_750"   , 4005, doSkim, doPrune );
  convertFiles(samples, "XZg_Spin0_ZToLL_W_0p014_M_900"   , 4006, doSkim, doPrune );
  convertFiles(samples, "XZg_Spin0_ZToLL_W_0p014_M_2000"  , 4008, doSkim, doPrune );

  std::cout << "-> Find your trees in: " << samples << "_conv" << std::endl;

  return 0;

}


void convertFiles( const std::string& samples, const std::string& expr, int id, bool skim, bool prune ) {

  std::string samplesFile = "../samples/samples_" + samples + ".dat";
  //std::string samplesFile = "../samples/samples_Run2015_25nsGolden_fromSnT.dat";

  std::string outdir = samples + "_conv";
  system( Form("mkdir -p %s", outdir.c_str()) );

  TChain* tree  = new TChain("mt2"); // this has all branches active: will get HLT/Flag info from this one
  TChain* tree2 = new TChain("mt2"); // this one has the branches deactivated: this will be cloned (so as to define new float branches)

  ifstream ifs(samplesFile);
  char buffer[500];
  int addedSamples = 0;

  while( ifs.getline(buffer, 500, '\n') ) {

    if (buffer[0] == '#') {
      continue; // Skip lines commented with '#'                                                                                                                                                                                 
    }

    std::string sampleFilePath(buffer);
    TString sampleFilePath_tstr(sampleFilePath);
    if( !(sampleFilePath_tstr.EndsWith(".root")) ) continue;

    if( !(sampleFilePath_tstr.Contains(expr)) ) continue;

    tree ->Add(sampleFilePath.c_str());
    tree2->Add(sampleFilePath.c_str());
    std::cout << "-> Added " << sampleFilePath << std::endl;
    addedSamples++;

  }

  if( addedSamples==0 ) {
    std::cout << "-> Didn't find any sample matching: " << expr << "... Skipping" << std::endl;
    return;
  }

  tree ->SetBranchStatus("evt_id",0);
  tree2->SetBranchStatus("evt_id",0);
  if( prune ) {
    tree->SetBranchStatus("genpart_*", 0);  
    tree->SetBranchStatus("genPart_*", 0);  
    tree->SetBranchStatus("genTau_*", 0);  
  }

  Int_t         snt_HLT_PFMET170;
  Int_t         snt_HLT_Photon90_R9Id90_HE10_IsoM;
  Int_t         snt_HLT_Photon75_R9Id90_HE10_IsoM;
  Int_t         snt_HLT_Photon120;
  Int_t         snt_HLT_Photon75;
  Int_t         snt_HLT_Photon165_HE10;
  Int_t         snt_HLT_Photon120_R9Id90_HE10_IsoM;
  Int_t         snt_HLT_Photon90;
  Int_t         snt_HLT_PFHT350_PFMET100;
  Int_t         snt_HLT_PFMETNoMu90_PFMHTNoMu90;
  Int_t         snt_HLT_PFMET90_PFMHT90;
  Int_t         snt_HLT_PFHT475_Prescale;
  Int_t         snt_HLT_PFHT350_Prescale;
  Int_t         snt_HLT_SingleMu;
  Int_t         snt_HLT_MuX_Ele12;
  Int_t         snt_HLT_Mu8_EleX;
  Int_t         snt_HLT_SingleEl;
  Int_t         snt_HLT_PFHT800;
  Int_t         snt_HLT_Photon155;
  Int_t         snt_HLT_PFHT900;
  Int_t         snt_HLT_Photon175;
  Int_t         snt_HLT_DiJet;
  Int_t         snt_HLT_DoubleEl;
  Int_t         snt_HLT_DoubleEl33;
  Int_t         snt_HLT_DoubleMu;
  Int_t         snt_Flag_EcalDeadCellTriggerPrimitiveFilter;
  Int_t         snt_Flag_trkPOG_manystripclus53X;
  Int_t         snt_Flag_ecalLaserCorrFilter;
  Int_t         snt_Flag_trkPOG_toomanystripclus53X;
  Int_t         snt_Flag_hcalLaserEventFilter;
  Int_t         snt_Flag_trkPOG_logErrorTooManyClusters;
  Int_t         snt_Flag_trkPOGFilters;
  Int_t         snt_Flag_trackingFailureFilter;
  Int_t         snt_Flag_CSCTightHaloFilter;
  Int_t         snt_Flag_HBHENoiseFilter;
  Int_t         snt_Flag_HBHEIsoNoiseFilter;
  Int_t         snt_Flag_goodVertices;
  Int_t         snt_Flag_METFilters;
  Int_t         snt_Flag_eeBadScFilter;

  TBranch        *b_snt_HLT_PFMET170;   //!
  TBranch        *b_snt_HLT_Photon90_R9Id90_HE10_IsoM;   //!
  TBranch        *b_snt_HLT_Photon75_R9Id90_HE10_IsoM;   //!
  TBranch        *b_snt_HLT_Photon120;   //!
  TBranch        *b_snt_HLT_Photon75;   //!
  TBranch        *b_snt_HLT_Photon165_HE10;   //!
  TBranch        *b_snt_HLT_Photon120_R9Id90_HE10_IsoM;   //!
  TBranch        *b_snt_HLT_Photon90;   //!
  TBranch        *b_snt_HLT_PFHT350_PFMET100;   //!
  TBranch        *b_snt_HLT_PFMETNoMu90_PFMHTNoMu90;   //!
  TBranch        *b_snt_HLT_PFMET90_PFMHT90;   //!
  TBranch        *b_snt_HLT_PFHT475_Prescale;   //!
  TBranch        *b_snt_HLT_PFHT350_Prescale;   //!
  TBranch        *b_snt_HLT_SingleMu;   //!
  TBranch        *b_snt_HLT_SingleEl;   //!
  TBranch        *b_snt_HLT_MuX_Ele12;   //!
  TBranch        *b_snt_HLT_Mu8_EleX;   //!
  TBranch        *b_snt_HLT_PFHT800;   //!
  TBranch        *b_snt_HLT_Photon155;   //!
  TBranch        *b_snt_HLT_PFHT900;   //!
  TBranch        *b_snt_HLT_Photon175;   //!
  TBranch        *b_snt_HLT_DiJet;   //!
  TBranch        *b_snt_HLT_DoubleEl;   //!
  TBranch        *b_snt_HLT_DoubleEl33;   //!
  TBranch        *b_snt_HLT_DoubleMu;   //!
  TBranch        *b_snt_Flag_EcalDeadCellTriggerPrimitiveFilter;   //!
  TBranch        *b_snt_Flag_trkPOG_manystripclus53X;   //!
  TBranch        *b_snt_Flag_ecalLaserCorrFilter;   //!
  TBranch        *b_snt_Flag_trkPOG_toomanystripclus53X;   //!
  TBranch        *b_snt_Flag_hcalLaserEventFilter;   //!
  TBranch        *b_snt_Flag_trkPOG_logErrorTooManyClusters;   //!
  TBranch        *b_snt_Flag_trkPOGFilters;   //!
  TBranch        *b_snt_Flag_trackingFailureFilter;   //!
  TBranch        *b_snt_Flag_CSCTightHaloFilter;   //!
  TBranch        *b_snt_Flag_HBHENoiseFilter;   //!
  TBranch        *b_snt_Flag_HBHEIsoNoiseFilter;   //!
  TBranch        *b_snt_Flag_goodVertices;   //!
  TBranch        *b_snt_Flag_METFilters;   //!
  TBranch        *b_snt_Flag_eeBadScFilter;   //!

  tree->SetBranchAddress("HLT_PFMET170", &snt_HLT_PFMET170, &b_snt_HLT_PFMET170);
  tree->SetBranchAddress("HLT_Photon90_R9Id90_HE10_IsoM", &snt_HLT_Photon90_R9Id90_HE10_IsoM, &b_snt_HLT_Photon90_R9Id90_HE10_IsoM);
  tree->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM", &snt_HLT_Photon75_R9Id90_HE10_IsoM, &b_snt_HLT_Photon75_R9Id90_HE10_IsoM);
  tree->SetBranchAddress("HLT_Photon120", &snt_HLT_Photon120, &b_snt_HLT_Photon120);
  tree->SetBranchAddress("HLT_Photon75", &snt_HLT_Photon75, &b_snt_HLT_Photon75);
  tree->SetBranchAddress("HLT_Photon165_HE10", &snt_HLT_Photon165_HE10, &b_snt_HLT_Photon165_HE10);
  tree->SetBranchAddress("HLT_Photon120_R9Id90_HE10_IsoM", &snt_HLT_Photon120_R9Id90_HE10_IsoM, &b_snt_HLT_Photon120_R9Id90_HE10_IsoM);
  tree->SetBranchAddress("HLT_Photon90", &snt_HLT_Photon90, &b_snt_HLT_Photon90);
  tree->SetBranchAddress("HLT_PFHT350_PFMET100", &snt_HLT_PFHT350_PFMET100, &b_snt_HLT_PFHT350_PFMET100);
  tree->SetBranchAddress("HLT_PFMETNoMu90_PFMHTNoMu90", &snt_HLT_PFMETNoMu90_PFMHTNoMu90, &b_snt_HLT_PFMETNoMu90_PFMHTNoMu90);
  tree->SetBranchAddress("HLT_PFMET90_PFMHT90", &snt_HLT_PFMET90_PFMHT90, &b_snt_HLT_PFMET90_PFMHT90);
  tree->SetBranchAddress("HLT_PFHT475_Prescale", &snt_HLT_PFHT475_Prescale, &b_snt_HLT_PFHT475_Prescale);
  tree->SetBranchAddress("HLT_PFHT350_Prescale", &snt_HLT_PFHT350_Prescale, &b_snt_HLT_PFHT350_Prescale);
  tree->SetBranchAddress("HLT_SingleMu", &snt_HLT_SingleMu, &b_snt_HLT_SingleMu);
  tree->SetBranchAddress("HLT_SingleEl", &snt_HLT_SingleEl, &b_snt_HLT_SingleEl);
  tree->SetBranchAddress("HLT_Mu8_EleX", &snt_HLT_Mu8_EleX, &b_snt_HLT_Mu8_EleX);
  tree->SetBranchAddress("HLT_MuX_Ele12", &snt_HLT_MuX_Ele12, &b_snt_HLT_MuX_Ele12);
  tree->SetBranchAddress("HLT_PFHT800", &snt_HLT_PFHT800, &b_snt_HLT_PFHT800);
  tree->SetBranchAddress("HLT_Photon155", &snt_HLT_Photon155, &b_snt_HLT_Photon155);
  tree->SetBranchAddress("HLT_PFHT900", &snt_HLT_PFHT900, &b_snt_HLT_PFHT900);
  tree->SetBranchAddress("HLT_Photon175", &snt_HLT_Photon175, &b_snt_HLT_Photon175);
  tree->SetBranchAddress("HLT_DiJet", &snt_HLT_DiJet, &b_snt_HLT_DiJet);
  tree->SetBranchAddress("HLT_DoubleEl", &snt_HLT_DoubleEl, &b_snt_HLT_DoubleEl);
  if( id<14000 )
    tree->SetBranchAddress("HLT_DoubleEl33", &snt_HLT_DoubleEl33, &b_snt_HLT_DoubleEl33);
  else
    tree->SetBranchAddress("HLT_DoubleEle33", &snt_HLT_DoubleEl33, &b_snt_HLT_DoubleEl33);
  tree->SetBranchAddress("HLT_DoubleMu", &snt_HLT_DoubleMu, &b_snt_HLT_DoubleMu);
  tree->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &snt_Flag_EcalDeadCellTriggerPrimitiveFilter, &b_snt_Flag_EcalDeadCellTriggerPrimitiveFilter);
  tree->SetBranchAddress("Flag_trkPOG_manystripclus53X", &snt_Flag_trkPOG_manystripclus53X, &b_snt_Flag_trkPOG_manystripclus53X);
  tree->SetBranchAddress("Flag_ecalLaserCorrFilter", &snt_Flag_ecalLaserCorrFilter, &b_snt_Flag_ecalLaserCorrFilter);
  tree->SetBranchAddress("Flag_trkPOG_toomanystripclus53X", &snt_Flag_trkPOG_toomanystripclus53X, &b_snt_Flag_trkPOG_toomanystripclus53X);
  tree->SetBranchAddress("Flag_hcalLaserEventFilter", &snt_Flag_hcalLaserEventFilter, &b_snt_Flag_hcalLaserEventFilter);
  tree->SetBranchAddress("Flag_trkPOG_logErrorTooManyClusters", &snt_Flag_trkPOG_logErrorTooManyClusters, &b_snt_Flag_trkPOG_logErrorTooManyClusters);
  tree->SetBranchAddress("Flag_trkPOGFilters", &snt_Flag_trkPOGFilters, &b_snt_Flag_trkPOGFilters);
  tree->SetBranchAddress("Flag_trackingFailureFilter", &snt_Flag_trackingFailureFilter, &b_snt_Flag_trackingFailureFilter);
  tree->SetBranchAddress("Flag_CSCTightHaloFilter", &snt_Flag_CSCTightHaloFilter, &b_snt_Flag_CSCTightHaloFilter);
  tree->SetBranchAddress("Flag_HBHENoiseFilter", &snt_Flag_HBHENoiseFilter, &b_snt_Flag_HBHENoiseFilter);
  tree->SetBranchAddress("Flag_HBHEIsoNoiseFilter", &snt_Flag_HBHEIsoNoiseFilter, &b_snt_Flag_HBHEIsoNoiseFilter);
  tree->SetBranchAddress("Flag_goodVertices", &snt_Flag_goodVertices, &b_snt_Flag_goodVertices);
  tree->SetBranchAddress("Flag_METFilters", &snt_Flag_METFilters, &b_snt_Flag_METFilters);
  tree->SetBranchAddress("Flag_eeBadScFilter", &snt_Flag_eeBadScFilter, &b_snt_Flag_eeBadScFilter);

  tree2->SetBranchStatus("HLT_PFMET170", 0 );
  tree2->SetBranchStatus("HLT_Photon90_R9Id90_HE10_IsoM", 0 );
  tree2->SetBranchStatus("HLT_Photon75_R9Id90_HE10_IsoM", 0 );
  tree2->SetBranchStatus("HLT_Photon120", 0 );
  tree2->SetBranchStatus("HLT_Photon75", 0 );
  tree2->SetBranchStatus("HLT_Photon165_HE10", 0 );
  tree2->SetBranchStatus("HLT_Photon120_R9Id90_HE10_IsoM", 0 );
  tree2->SetBranchStatus("HLT_Photon90", 0 );
  tree2->SetBranchStatus("HLT_PFHT350_PFMET100", 0 );
  tree2->SetBranchStatus("HLT_PFMETNoMu90_PFMHTNoMu90", 0 );
  tree2->SetBranchStatus("HLT_PFMET90_PFMHT90", 0 );
  tree2->SetBranchStatus("HLT_PFHT475_Prescale", 0 );
  tree2->SetBranchStatus("HLT_PFHT350_Prescale", 0 );
  tree2->SetBranchStatus("HLT_SingleMu", 0 );
  tree2->SetBranchStatus("HLT_SingleEl", 0 );
  tree2->SetBranchStatus("HLT_Mu8_EleX", 0 );
  tree2->SetBranchStatus("HLT_MuX_Ele12", 0 );
  tree2->SetBranchStatus("HLT_PFHT800", 0 );
  tree2->SetBranchStatus("HLT_Photon155", 0 );
  tree2->SetBranchStatus("HLT_PFHT900", 0 );
  tree2->SetBranchStatus("HLT_Photon175", 0 );
  tree2->SetBranchStatus("HLT_DiJet", 0 );
  tree2->SetBranchStatus("HLT_DoubleEl", 0 );
  if( id<14000 )
    tree2->SetBranchStatus("HLT_DoubleEl33", 0 );
  else
    tree2->SetBranchStatus("HLT_DoubleEle33", 0 );
  tree2->SetBranchStatus("HLT_DoubleMu", 0 );
  tree2->SetBranchStatus("Flag_EcalDeadCellTriggerPrimitiveFilter", 0 );
  tree2->SetBranchStatus("Flag_trkPOG_manystripclus53X", 0 );
  tree2->SetBranchStatus("Flag_ecalLaserCorrFilter", 0 );
  tree2->SetBranchStatus("Flag_trkPOG_toomanystripclus53X", 0 );
  tree2->SetBranchStatus("Flag_hcalLaserEventFilter", 0 );
  tree2->SetBranchStatus("Flag_trkPOG_logErrorTooManyClusters", 0 );
  tree2->SetBranchStatus("Flag_trkPOGFilters", 0 );
  tree2->SetBranchStatus("Flag_trackingFailureFilter", 0 );
  tree2->SetBranchStatus("Flag_CSCTightHaloFilter", 0 );
  tree2->SetBranchStatus("Flag_HBHENoiseFilter", 0 );
  tree2->SetBranchStatus("Flag_HBHEIsoNoiseFilter", 0 );
  tree2->SetBranchStatus("Flag_goodVertices", 0 );
  tree2->SetBranchStatus("Flag_METFilters", 0 );
  tree2->SetBranchStatus("Flag_eeBadScFilter", 0 );



  std::string outFileName = outdir + "/" + expr;
  if( skim ) outFileName = outFileName + "_skim";
  if( loose_skim ) outFileName = outFileName + "_loose";
  outFileName = outFileName + ".root";
  TFile* outFile = TFile::Open( outFileName.c_str(), "recreate" );
  outFile->cd();

  TTree* outTree = tree2->CloneTree(0);


  Float_t         HLT_PFMET170;
  Float_t         HLT_Photon90_R9Id90_HE10_IsoM;
  Float_t         HLT_Photon75_R9Id90_HE10_IsoM;
  Float_t         HLT_Photon120;
  Float_t         HLT_Photon75;
  Float_t         HLT_Photon165_HE10;
  Float_t         HLT_Photon120_R9Id90_HE10_IsoM;
  Float_t         HLT_Photon90;
  Float_t         HLT_PFHT350_PFMET100;
  Float_t         HLT_PFMETNoMu90_PFMHTNoMu90;
  Float_t         HLT_PFMET90_PFMHT90;
  Float_t         HLT_PFHT475_Prescale;
  Float_t         HLT_PFHT350_Prescale;
  Float_t         HLT_SingleMu;
  Float_t         HLT_MuX_Ele12;
  Float_t         HLT_Mu8_EleX;
  Float_t         HLT_SingleEl;
  Float_t         HLT_PFHT800;
  Float_t         HLT_Photon155;
  Float_t         HLT_PFHT900;
  Float_t         HLT_Photon175;
  Float_t         HLT_DiJet;
  Float_t         HLT_DoubleEl;
  Float_t         HLT_DoubleEl33;
  Float_t         HLT_DoubleMu;
  Float_t         Flag_EcalDeadCellTriggerPrimitiveFilter;
  Float_t         Flag_trkPOG_manystripclus53X;
  Float_t         Flag_ecalLaserCorrFilter;
  Float_t         Flag_trkPOG_toomanystripclus53X;
  Float_t         Flag_hcalLaserEventFilter;
  Float_t         Flag_trkPOG_logErrorTooManyClusters;
  Float_t         Flag_trkPOGFilters;
  Float_t         Flag_trackingFailureFilter;
  Float_t         Flag_CSCTightHaloFilter;
  Float_t         Flag_HBHENoiseFilter;
  Float_t         Flag_HBHEIsoNoiseFilter;
  Float_t         Flag_goodVertices;
  Float_t         Flag_METFilters;
  Float_t         Flag_eeBadScFilter;

  outTree->Branch("HLT_PFMET170", &HLT_PFMET170, "HLT_PFMET170/F" );
  outTree->Branch("HLT_Photon90_R9Id90_HE10_IsoM", &HLT_Photon90_R9Id90_HE10_IsoM, "HLT_Photon90_R9Id90_HE10_IsoM/F" );
  outTree->Branch("HLT_Photon75_R9Id90_HE10_IsoM", &HLT_Photon75_R9Id90_HE10_IsoM, "HLT_Photon75_R9Id90_HE10_IsoM/F" );
  outTree->Branch("HLT_Photon120", &HLT_Photon120, "HLT_Photon120/F" );
  outTree->Branch("HLT_Photon75", &HLT_Photon75, "HLT_Photon75/F" );
  outTree->Branch("HLT_Photon165_HE10", &HLT_Photon165_HE10, "HLT_Photon165_HE10/F" );
  outTree->Branch("HLT_Photon120_R9Id90_HE10_IsoM", &HLT_Photon120_R9Id90_HE10_IsoM, "HLT_Photon120_R9Id90_HE10_IsoM/F" );
  outTree->Branch("HLT_Photon90", &HLT_Photon90, "HLT_Photon90/F" );
  outTree->Branch("HLT_PFHT350_PFMET100", &HLT_PFHT350_PFMET100, "HLT_PFHT350_PFMET100/F" );
  outTree->Branch("HLT_PFMETNoMu90_PFMHTNoMu90", &HLT_PFMETNoMu90_PFMHTNoMu90, "HLT_PFMETNoMu90_PFMHTNoMu90/F" );
  outTree->Branch("HLT_PFMET90_PFMHT90", &HLT_PFMET90_PFMHT90, "HLT_PFMET90_PFMHT90/F" );
  outTree->Branch("HLT_PFHT475_Prescale", &HLT_PFHT475_Prescale, "HLT_PFHT475_Prescale/F" );
  outTree->Branch("HLT_PFHT350_Prescale", &HLT_PFHT350_Prescale, "HLT_PFHT350_Prescale/F" );
  outTree->Branch("HLT_SingleMu", &HLT_SingleMu, "HLT_SingleMu/F" );
  outTree->Branch("HLT_SingleEl", &HLT_SingleEl, "HLT_SingleEl/F" );
  outTree->Branch("HLT_Mu8_EleX", &HLT_Mu8_EleX, "HLT_Mu8_EleX/F" );
  outTree->Branch("HLT_MuX_Ele12", &HLT_MuX_Ele12, "HLT_MuX_Ele12/F" );
  outTree->Branch("HLT_PFHT800", &HLT_PFHT800, "HLT_PFHT800/F" );
  outTree->Branch("HLT_Photon155", &HLT_Photon155, "HLT_Photon155/F" );
  outTree->Branch("HLT_PFHT900", &HLT_PFHT900, "HLT_PFHT900/F" );
  outTree->Branch("HLT_Photon175", &HLT_Photon175, "HLT_Photon175/F" );
  outTree->Branch("HLT_DiJet", &HLT_DiJet, "HLT_DiJet/F" );
  outTree->Branch("HLT_DoubleEl", &HLT_DoubleEl, "HLT_DoubleEl/F" );
  outTree->Branch("HLT_DoubleEl33", &HLT_DoubleEl33, "HLT_DoubleEl33/F" );
  outTree->Branch("HLT_DoubleMu", &HLT_DoubleMu, "HLT_DoubleMu/F" );
  outTree->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, "Flag_EcalDeadCellTriggerPrimitiveFilter/F" );
  outTree->Branch("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X, "Flag_trkPOG_manystripclus53X/F" );
  outTree->Branch("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter, "Flag_ecalLaserCorrFilter/F" );
  outTree->Branch("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X, "Flag_trkPOG_toomanystripclus53X/F" );
  outTree->Branch("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter, "Flag_hcalLaserEventFilter/F" );
  outTree->Branch("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters, "Flag_trkPOG_logErrorTooManyClusters/F" );
  outTree->Branch("Flag_trkPOGFilters", &Flag_trkPOGFilters, "Flag_trkPOGFilters/F" );
  outTree->Branch("Flag_trackingFailureFilter", &Flag_trackingFailureFilter, "Flag_trackingFailureFilter/F" );
  outTree->Branch("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter, "Flag_CSCTightHaloFilter/F" );
  outTree->Branch("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, "Flag_HBHENoiseFilter/F" );
  outTree->Branch("Flag_HBHEIsoNoiseFilter", &Flag_HBHEIsoNoiseFilter, "Flag_HBHEIsoNoiseFilter/F" );
  outTree->Branch("Flag_goodVertices", &Flag_goodVertices, "Flag_goodVertices/F" );
  outTree->Branch("Flag_METFilters", &Flag_METFilters, "Flag_METFilters/F" );
  outTree->Branch("Flag_eeBadScFilter", &Flag_eeBadScFilter, "Flag_eeBadScFilter/F" );

  int id_branch;
  outTree->Branch("evt_id", &id_branch );


  // for the skimming:
  Int_t           ngamma;
  Float_t         gamma_pt[20];   //[ngamma]
  Float_t         gamma_eta[20];   //[ngamma]
  Float_t         gamma_phi[20];   //[ngamma]
  Float_t         gamma_mass[20];   //[ngamma]
  Float_t         gamma_chHadIso[20];   //[ngamma]
  TBranch        *b_ngamma;   //!
  TBranch        *b_gamma_pt;   //!
  TBranch        *b_gamma_eta;   //!
  TBranch        *b_gamma_phi;   //!
  TBranch        *b_gamma_mass;   //!
  TBranch        *b_gamma_chHadIso;   //!
  tree->SetBranchAddress("ngamma", &ngamma, &b_ngamma);
  tree->SetBranchAddress("gamma_pt", gamma_pt, &b_gamma_pt);
  tree->SetBranchAddress("gamma_eta", gamma_eta, &b_gamma_eta);
  tree->SetBranchAddress("gamma_phi", gamma_phi, &b_gamma_phi);
  tree->SetBranchAddress("gamma_mass", gamma_mass, &b_gamma_mass);
  tree->SetBranchAddress("gamma_chHadIso", gamma_chHadIso, &b_gamma_chHadIso);

  Int_t           nlep;
  Float_t         lep_pt[20];   //[nlep]
  Float_t         lep_eta[20];   //[nlep]
  Float_t         lep_phi[20];   //[nlep]
  Float_t         lep_mass[20];   //[nlep]
  TBranch        *b_nlep;   //!
  TBranch        *b_lep_pt;   //!
  TBranch        *b_lep_eta;   //!
  TBranch        *b_lep_phi;   //!
  TBranch        *b_lep_mass;   //!
  tree->SetBranchAddress("nlep", &nlep, &b_nlep);
  tree->SetBranchAddress("lep_pt", lep_pt, &b_lep_pt);
  tree->SetBranchAddress("lep_eta", lep_eta, &b_lep_eta);
  tree->SetBranchAddress("lep_phi", lep_phi, &b_lep_phi);
  tree->SetBranchAddress("lep_mass", lep_mass, &b_lep_mass);


  int nentries = tree->GetEntries();

  for( int iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 100000 == 0 ) std::cout << " Entry: " << iEntry << " / " << nentries << std::endl;

    tree ->GetEntry(iEntry);
    tree2->GetEntry(iEntry);

    id_branch = id;

    HLT_PFMET170                            = (float)snt_HLT_PFMET170;
    HLT_Photon90_R9Id90_HE10_IsoM           = (float)snt_HLT_Photon90_R9Id90_HE10_IsoM;
    HLT_Photon75_R9Id90_HE10_IsoM           = (float)snt_HLT_Photon75_R9Id90_HE10_IsoM;
    HLT_Photon120                           = (float)snt_HLT_Photon120;
    HLT_Photon75                            = (float)snt_HLT_Photon75;
    HLT_Photon165_HE10                      = (float)snt_HLT_Photon165_HE10;
    HLT_Photon120_R9Id90_HE10_IsoM          = (float)snt_HLT_Photon120_R9Id90_HE10_IsoM;
    HLT_Photon90                            = (float)snt_HLT_Photon90;
    HLT_PFHT350_PFMET100                    = (float)snt_HLT_PFHT350_PFMET100;
    HLT_PFMETNoMu90_PFMHTNoMu90             = (float)snt_HLT_PFMETNoMu90_PFMHTNoMu90;
    HLT_PFMET90_PFMHT90                     = (float)snt_HLT_PFMET90_PFMHT90;
    HLT_PFHT475_Prescale                    = (float)snt_HLT_PFHT475_Prescale;
    HLT_PFHT350_Prescale                    = (float)snt_HLT_PFHT350_Prescale;
    HLT_SingleMu                            = (float)snt_HLT_SingleMu;
    HLT_MuX_Ele12                           = (float)snt_HLT_MuX_Ele12;
    HLT_Mu8_EleX                            = (float)snt_HLT_Mu8_EleX;
    HLT_SingleEl                            = (float)snt_HLT_SingleEl;
    HLT_PFHT800                             = (float)snt_HLT_PFHT800;
    HLT_Photon155                           = (float)snt_HLT_Photon155;
    HLT_PFHT900                             = (float)snt_HLT_PFHT900;
    HLT_Photon175                           = (float)snt_HLT_Photon175;
    HLT_DiJet                               = (float)snt_HLT_DiJet;
    HLT_DoubleEl                            = (float)snt_HLT_DoubleEl;
    HLT_DoubleEl33                          = (float)snt_HLT_DoubleEl33;
    HLT_DoubleMu                            = (float)snt_HLT_DoubleMu;
    Flag_EcalDeadCellTriggerPrimitiveFilter = (float)snt_Flag_EcalDeadCellTriggerPrimitiveFilter;
    Flag_trkPOG_manystripclus53X            = (float)snt_Flag_trkPOG_manystripclus53X;
    Flag_ecalLaserCorrFilter                = (float)snt_Flag_ecalLaserCorrFilter;
    Flag_trkPOG_toomanystripclus53X         = (float)snt_Flag_trkPOG_toomanystripclus53X;
    Flag_hcalLaserEventFilter               = (float)snt_Flag_hcalLaserEventFilter;
    Flag_trkPOG_logErrorTooManyClusters     = (float)snt_Flag_trkPOG_logErrorTooManyClusters;
    Flag_trkPOGFilters                      = (float)snt_Flag_trkPOGFilters;
    Flag_trackingFailureFilter              = (float)snt_Flag_trackingFailureFilter;
    Flag_CSCTightHaloFilter                 = (float)snt_Flag_CSCTightHaloFilter;
    Flag_HBHENoiseFilter                    = (float)snt_Flag_HBHENoiseFilter;
    Flag_HBHEIsoNoiseFilter                 = (float)snt_Flag_HBHEIsoNoiseFilter;
    Flag_goodVertices                       = (float)snt_Flag_goodVertices;
    Flag_METFilters                         = (float)snt_Flag_METFilters;
    Flag_eeBadScFilter                      = (float)snt_Flag_eeBadScFilter;


    if( skim ) {

      if( !(Flag_HBHENoiseFilter>0. && Flag_HBHEIsoNoiseFilter>0. && Flag_eeBadScFilter>0.) ) continue;
      if( !loose_skim )
        if( ngamma==0 ) continue;
      if( nlep!=2 ) continue;
      if( lep_pt[0]<25. || lep_pt[1]<20. ) continue;

      //int foundPhot = -1;
      //for( int iph=0; iph<ngamma; ++iph ) {
      //  if( gamma_chHadIso[iph]>10. ) continue; // loose iso
      //  TLorentzVector phot;
      //  phot.SetPtEtaPhiM( gamma_pt[iph], gamma_eta[iph], gamma_phi[iph], gamma_mass[iph] );
      //  for( int il=0; il<nlep; ++il ) {
      //    TLorentzVector lep;
      //    lep.SetPtEtaPhiM( lep_pt[il], lep_eta[il], lep_phi[il], lep_mass[il] );
      //    if( phot.DeltaR(lep)<0.3 ) continue;
      //  }
      //  foundPhot = iph;
      //  break;
      //}

      //if( foundPhot<0 ) continue;
      //if( gamma_pt[foundPhot]<40. ) continue;

    } // skim
         
    outTree->Fill();

  }

  outTree->Write();

  outFile->Close();

}
