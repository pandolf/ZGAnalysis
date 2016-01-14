#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"


void convertFiles( const std::string& samplesFile, const std::string& expr, int id, bool skim=false, bool prune=false, const std::string& outdir="." );



int main() {


  std::string samplesFile = "../samples/samples_Run2015_25nsGolden_fromSnT.dat";

  bool doSkim = true;
  bool doPrune = true;

  convertFiles(samplesFile, "DoubleEG"    , 4, doSkim, doPrune);
  convertFiles(samplesFile, "DoubleMuon"  , 5, doSkim, doPrune);
  convertFiles(samplesFile, "SinglePhoton", 7, doSkim, doPrune);

  return 0;

}


void convertFiles( const std::string& samplesFile, const std::string& expr, int id, bool skim, bool prune, const std::string& outdir ) {

  system( Form("mkdir -p %s", outdir.c_str()) );

  TChain* tree = new TChain("mt2");

  ifstream ifs(samplesFile);
  char buffer[500];

  while( ifs.getline(buffer, 500, '\n') ) {

    if (buffer[0] == '#') {
      continue; // Skip lines commented with '#'                                                                                                                                                                                 
    }

    std::string sampleFilePath(buffer);
    TString sampleFilePath_tstr(sampleFilePath);
    if( !(sampleFilePath_tstr.EndsWith(".root")) ) continue;

    if( !(sampleFilePath_tstr.Contains(expr)) ) continue;

    tree->Add(sampleFilePath.c_str());
    std::cout << "-> Added " << sampleFilePath << std::endl;

  }

  // disable Flag and HLT as i want to convert to float
  tree->SetBranchStatus("Flag_*",0);
  tree->SetBranchStatus("HLT_*",0);
  tree->SetBranchStatus("evt_id",0);
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



  TFile* outFile = TFile::Open( Form("%s/%s.root", outdir.c_str(), expr.c_str()), "recreate" );
  outFile->cd();

  TTree* outTree = tree->CloneTree(0);


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
  TBranch        *b_ngamma;   //!
  TBranch        *b_gamma_pt;   //!
  tree->SetBranchAddress("ngamma", &ngamma, &b_ngamma);
  tree->SetBranchAddress("gamma_pt", gamma_pt, &b_gamma_pt);


  int nentries = tree->GetEntries();

  for( int iEntry=0; iEntry<nentries; ++iEntry ) {

    tree->GetEntry(iEntry);

    id_branch = id;

    HLT_PFMET170                            = snt_HLT_PFMET170;
    HLT_Photon90_R9Id90_HE10_IsoM           = snt_HLT_Photon90_R9Id90_HE10_IsoM;
    HLT_Photon75_R9Id90_HE10_IsoM           = snt_HLT_Photon75_R9Id90_HE10_IsoM;
    HLT_Photon120                           = snt_HLT_Photon120;
    HLT_Photon75                            = snt_HLT_Photon75;
    HLT_Photon165_HE10                      = snt_HLT_Photon165_HE10;
    HLT_Photon120_R9Id90_HE10_IsoM          = snt_HLT_Photon120_R9Id90_HE10_IsoM;
    HLT_Photon90                            = snt_HLT_Photon90;
    HLT_PFHT350_PFMET100                    = snt_HLT_PFHT350_PFMET100;
    HLT_PFMETNoMu90_PFMHTNoMu90             = snt_HLT_PFMETNoMu90_PFMHTNoMu90;
    HLT_PFMET90_PFMHT90                     = snt_HLT_PFMET90_PFMHT90;
    HLT_PFHT475_Prescale                    = snt_HLT_PFHT475_Prescale;
    HLT_PFHT350_Prescale                    = snt_HLT_PFHT350_Prescale;
    HLT_SingleMu                            = snt_HLT_SingleMu;
    HLT_MuX_Ele12                           = snt_HLT_MuX_Ele12;
    HLT_Mu8_EleX                            = snt_HLT_Mu8_EleX;
    HLT_SingleEl                            = snt_HLT_SingleEl;
    HLT_PFHT800                             = snt_HLT_PFHT800;
    HLT_Photon155                           = snt_HLT_Photon155;
    HLT_PFHT900                             = snt_HLT_PFHT900;
    HLT_Photon175                           = snt_HLT_Photon175;
    HLT_DiJet                               = snt_HLT_DiJet;
    HLT_DoubleEl                            = snt_HLT_DoubleEl;
    HLT_DoubleMu                            = snt_HLT_DoubleMu;
    Flag_EcalDeadCellTriggerPrimitiveFilter = snt_Flag_EcalDeadCellTriggerPrimitiveFilter;
    Flag_trkPOG_manystripclus53X            = snt_Flag_trkPOG_manystripclus53X;
    Flag_ecalLaserCorrFilter                = snt_Flag_ecalLaserCorrFilter;
    Flag_trkPOG_toomanystripclus53X         = snt_Flag_trkPOG_toomanystripclus53X;
    Flag_hcalLaserEventFilter               = snt_Flag_hcalLaserEventFilter;
    Flag_trkPOG_logErrorTooManyClusters     = snt_Flag_trkPOG_logErrorTooManyClusters;
    Flag_trkPOGFilters                      = snt_Flag_trkPOGFilters;
    Flag_trackingFailureFilter              = snt_Flag_trackingFailureFilter;
    Flag_CSCTightHaloFilter                 = snt_Flag_CSCTightHaloFilter;
    Flag_HBHENoiseFilter                    = snt_Flag_HBHENoiseFilter;
    Flag_HBHEIsoNoiseFilter                 = snt_Flag_HBHEIsoNoiseFilter;
    Flag_goodVertices                       = snt_Flag_goodVertices;
    Flag_METFilters                         = snt_Flag_METFilters;
    Flag_eeBadScFilter                      = snt_Flag_eeBadScFilter;

    if( skim ) 
      if( !(ngamma>0 && gamma_pt[0]>30. ) ) continue;
         
    outTree->Fill();

  }

  outTree->Write();

  outFile->Close();

}
