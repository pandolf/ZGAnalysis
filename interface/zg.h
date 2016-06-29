//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Dec 12 23:39:08 2014 by ROOT version 5.34/10
// from TTree mt2/treeProducerSusyFullHad
// found on file: /scratch/mmasciov/12December/POSTPROCESS/TTJets_MSDecaysCKM_central_PU_S14_POSTLS170_babytree_post.root
//////////////////////////////////////////////////////////

#ifndef zg_h
#define zg_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class ZGTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           lumi;
   ULong64_t       evt;
   Int_t           isData;
   Int_t           isGolden;
   Int_t           isSilver;
   Float_t         xsec;
   Int_t         HLT_PFMET170;
   Int_t         HLT_Photon90_R9Id90_HE10_IsoM;
   Int_t         HLT_Photon75_R9Id90_HE10_IsoM;
   Int_t         HLT_Photon120;
   Int_t         HLT_Photon75;
   Int_t         HLT_Photon165_HE10;
   Int_t         HLT_Photon120_R9Id90_HE10_IsoM;
   Int_t         HLT_Photon90;
   Int_t         HLT_PFHT350_PFMET100;
   Int_t         HLT_PFMETNoMu90_PFMHTNoMu90;
   Int_t         HLT_PFMET90_PFMHT90;
   Int_t         HLT_PFHT475_Prescale;
   Int_t         HLT_PFHT350_Prescale;
   Int_t         HLT_SingleMu;
   Int_t         HLT_MuX_Ele12;
   Int_t         HLT_Mu8_EleX;
   Int_t         HLT_SingleEl;
   Int_t         HLT_PFHT800;
   Int_t         HLT_Photon155;
   Int_t         HLT_PFHT900;
   Int_t         HLT_Photon175;
   Int_t         HLT_DoubleEl;
   Int_t         HLT_DoubleEle33;
   Int_t         HLT_DoubleMu;
   Int_t         HLT_Mu30_TkMu11;
   //Float_t         HLT_PFMET170;
   //Float_t         HLT_Photon90_R9Id90_HE10_IsoM;
   //Float_t         HLT_Photon75_R9Id90_HE10_IsoM;
   //Float_t         HLT_Photon120;
   //Float_t         HLT_Photon75;
   //Float_t         HLT_Photon165_HE10;
   //Float_t         HLT_Photon120_R9Id90_HE10_IsoM;
   //Float_t         HLT_Photon90;
   //Float_t         HLT_PFHT350_PFMET100;
   //Float_t         HLT_PFMETNoMu90_PFMHTNoMu90;
   //Float_t         HLT_PFMET90_PFMHT90;
   //Float_t         HLT_PFHT475_Prescale;
   //Float_t         HLT_PFHT350_Prescale;
   //Float_t         HLT_SingleMu;
   //Float_t         HLT_MuX_Ele12;
   //Float_t         HLT_Mu8_EleX;
   //Float_t         HLT_SingleEl;
   //Float_t         HLT_PFHT800;
   //Float_t         HLT_Photon155;
   //Float_t         HLT_PFHT900;
   //Float_t         HLT_Photon175;
   //Float_t         HLT_DoubleEl;
   //Float_t         HLT_DoubleEle33;
   //Float_t         HLT_DoubleMu;
   //Float_t         HLT_Mu30_TkMu11;
   //Int_t           Flag_EcalDeadCellTriggerPrimitiveFilter;
   //Int_t           Flag_trkPOG_manystripclus53X;
   //Int_t           Flag_ecalLaserCorrFilter;
   //Int_t           Flag_trkPOG_toomanystripclus53X;
   //Int_t           Flag_hcalLaserEventFilter;
   //Int_t           Flag_trkPOG_logErrorTooManyClusters;
   //Int_t           Flag_trkPOGFilters;
   //Int_t           Flag_trackingFailureFilter;
   //Int_t           Flag_CSCTightHaloFilter;
   //Int_t           Flag_HBHENoiseFilter;
   //Int_t           Flag_HBHENoiseIsoFilter;
   //Int_t           Flag_goodVertices;
   //Int_t           Flag_METFilters;
   //Int_t           Flag_eeBadScFilter;
   Int_t           Flag_EcalDeadCellTriggerPrimitiveFilter;
   Int_t           Flag_HBHENoiseIsoFilter;
   Int_t           Flag_goodVertices;
   Int_t           Flag_HBHENoiseFilter;
   Int_t           Flag_CSCTightHalo2015Filter;
   Int_t           Flag_eeBadScFilter;
   Int_t           Flag_globalTightHalo2016Filter;
   Float_t         puWeight;
   Float_t         nTrueInt;
   Float_t         genWeight;
   Float_t         rho;
   Int_t           nVert;
   Float_t         genQScale;
   Float_t         met_rawPt;
   Float_t         met_rawPhi;
   Int_t           nMuons10;
   Int_t           nElectrons10;
   Int_t           nTaus20;
   Int_t           nGammas20;
   Int_t           nLepLowMT;
   Int_t           nPFLep5LowMT;
   Int_t           nPFHad10LowMT;
   Float_t         met_pt;
   Float_t         met_eta;
   Float_t         met_phi;
   Float_t         met_mass;
   Float_t         met_sumEt;
   Float_t         met_genPt;
   Float_t         met_genPhi;
   Float_t         met_genEta;
   Int_t           nGenPart;
   Int_t           GenPart_motherId[11];   //[nGenPart]
   Int_t           GenPart_grandmotherId[11];   //[nGenPart]
   Int_t           GenPart_sourceId[11];   //[nGenPart]
   Float_t         GenPart_charge[11];   //[nGenPart]
   Int_t           GenPart_status[11];   //[nGenPart]
   Int_t           GenPart_pdgId[11];   //[nGenPart]
   Float_t         GenPart_pt[11];   //[nGenPart]
   Float_t         GenPart_eta[11];   //[nGenPart]
   Float_t         GenPart_phi[11];   //[nGenPart]
   Float_t         GenPart_mass[11];   //[nGenPart]
   Int_t           GenPart_motherIndex[11];   //[nGenPart]
   Int_t           ngenPart;
   Float_t         genPart_pt[200];   //[ngenPart]
   Float_t         genPart_eta[200];   //[ngenPart]
   Float_t         genPart_phi[200];   //[ngenPart]
   Float_t         genPart_mass[200];   //[ngenPart]
   Int_t           genPart_pdgId[200];   //[ngenPart]
   Float_t         genPart_charge[200];   //[ngenPart]
   Int_t           genPart_status[200];   //[ngenPart]
   Int_t           genPart_motherId[200];   //[ngenPart]
   Int_t           genPart_motherIndex[200];   //[ngenPart]
   Int_t           genPart_grandmaId[200];   //[ngenPart]
   Int_t           njet;
   Float_t         jet_pt[50];   //[njet]
   Float_t         jet_eta[50];   //[njet]
   Float_t         jet_phi[50];   //[njet]
   Float_t         jet_mass[50];   //[njet]
   Float_t         jet_btagCSV[50];   //[njet]
   Float_t         jet_rawPt[50];   //[njet]
   Float_t         jet_mcPt[50];   //[njet]
   Int_t           jet_mcFlavour[50];   //[njet]
   Int_t           jet_mcMatchId[50];   //[njet]
   Int_t           jet_mcMatchFlav[50];   //[njet]
   Float_t         jet_area[50];   //[njet]
   Int_t           jet_puId[50];   //[njet]
   Int_t           jet_id[50];   //[njet]
   Float_t         jet_ptd[50];   //[njet]
   Float_t         jet_axis2[50];   //[njet]
   Int_t           jet_mult[50];   //[njet]
   Float_t         jet_qgl[50];   //[njet]
   Int_t           jet_partonId[50];   //[njet]
   Int_t           jet_partonMotherId[50];   //[njet]
   Float_t         jet_chHEF[50];   //[njet]
   Float_t         jet_neHEF[50];   //[njet]
   Float_t         jet_phEF[50];   //[njet]
   Float_t         jet_eEF[50];   //[njet]
   Float_t         jet_muEF[50];   //[njet]
   Float_t         jet_HFHEF[50];   //[njet]
   Float_t         jet_HFEMEF[50];   //[njet]
   Int_t           jet_chHMult[50];   //[njet]
   Int_t           jet_neHMult[50];   //[njet]
   Int_t           jet_phMult[50];   //[njet]
   Int_t           jet_eMult[50];   //[njet]
   Int_t           jet_muMult[50];   //[njet]
   Int_t           jet_HFHMult[50];   //[njet]
   Int_t           jet_HFEMMult[50];   //[njet]
   Int_t           nlep;
   Float_t         lep_pt[5];   //[nlep]
   Float_t         lep_eta[5];   //[nlep]
   Float_t         lep_phi[5];   //[nlep]
   Float_t         lep_mass[5];   //[nlep]
   Int_t           lep_pdgId[5];   //[nlep]
   Int_t           lep_charge[5];   //[nlep]
   Float_t         lep_dxy[5];   //[nlep]
   Float_t         lep_dz[5];   //[nlep]
   Float_t         lep_edxy[5];   //[nlep]
   Float_t         lep_edz[5];   //[nlep]
   Float_t         lep_ip3d[5];   //[nlep]
   Float_t         lep_sip3d[5];   //[nlep]
   Int_t           lep_tightId[5];   //[nlep]
   Int_t           lep_convVeto[5];   //[nlep]
   Int_t           lep_lostHits[5];   //[nlep]
   Int_t           lep_looseIdSusy[5];   //[nlep]
   Float_t         lep_relIso03[5];   //[nlep]
   Float_t         lep_relIso04[5];   //[nlep]
   Float_t         lep_miniRelIso[5];   //[nlep]
   Float_t         lep_chargedHadRelIso03[5];   //[nlep]
   Float_t         lep_chargedHadRelIso04[5];   //[nlep]
   Int_t           lep_convVetoFull[5];   //[nlep]
   Int_t           lep_eleCutIdCSA14_25ns_v1[5];   //[nlep]
   Int_t           lep_eleCutIdCSA14_50ns_v1[5];   //[nlep]
   Int_t           lep_eleMVAId[5];   //[nlep]
   Int_t           lep_tightCharge[5];   //[nlep]
   Float_t         lep_mvaId[5];   //[nlep]
   Float_t         lep_mvaIdTrig[5];   //[nlep]
   Float_t         lep_nStations[5];   //[nlep]
   Float_t         lep_trkKink[5];   //[nlep]
   Float_t         lep_caloCompatibility[5];   //[nlep]
   Float_t         lep_globalTrackChi2[5];   //[nlep]
   Int_t           lep_trackerLayers[5];   //[nlep]
   Int_t           lep_pixelLayers[5];   //[nlep]
   Float_t         lep_mvaTTH[5];   //[nlep]
   Float_t         lep_mvaSusy[5];   //[nlep]
   Float_t         lep_jetPtRatio[5];   //[nlep]
   Float_t         lep_jetBTagCSV[5];   //[nlep]
   Float_t         lep_jetDR[5];   //[nlep]
   Int_t           lep_mcMatchId[5];   //[nlep]
   Int_t           lep_mcMatchAny[5];   //[nlep]
   Int_t           lep_mcMatchTau[5];   //[nlep]
   Float_t         lep_r9[5];   //[nlep]
   Int_t           ngenLepFromTau;
   Float_t         genLepFromTau_pt[5];   //[ngenLepFromTau]
   Float_t         genLepFromTau_eta[5];   //[ngenLepFromTau]
   Float_t         genLepFromTau_phi[5];   //[ngenLepFromTau]
   Float_t         genLepFromTau_mass[5];   //[ngenLepFromTau]
   Int_t           genLepFromTau_pdgId[5];   //[ngenLepFromTau]
   Float_t         genLepFromTau_charge[5];   //[ngenLepFromTau]
   Int_t           genLepFromTau_status[5];   //[ngenLepFromTau]
   Int_t           genLepFromTau_sourceId[5];   //[ngenLepFromTau]
   Int_t           ngamma;
   Float_t         gamma_pt[20];   //[ngamma]
   Float_t         gamma_eta[20];   //[ngamma]
   Float_t         gamma_phi[20];   //[ngamma]
   Float_t         gamma_mass[20];   //[ngamma]
   Int_t           gamma_pdgId[20];   //[ngamma]
   Int_t           gamma_idCutBased[20];   //[ngamma]
   Float_t         gamma_hOverE[20];   //[ngamma]
   Float_t         gamma_r9[20];   //[ngamma]
   Float_t         gamma_sigmaIetaIeta[20];   //[ngamma]
   Float_t         gamma_chHadIso[20];   //[ngamma]
   Float_t         gamma_chHadIso04[20];   //[ngamma]
   Float_t         gamma_chHadIsoRC[20];   //[ngamma]
   Float_t         gamma_chHadIsoRC04[20];   //[ngamma]
   Float_t         gamma_neuHadIso[20];   //[ngamma]
   Float_t         gamma_phIso[20];   //[ngamma]
   Int_t           gamma_mcMatchId[20];   //[ngamma]
   Float_t         gamma_mcPt[20];   //[ngamma]
   Float_t         gamma_genIso03[20];   //[ngamma]
   Float_t         gamma_genIso04[20];   //[ngamma]
   Float_t         gamma_drMinParton[20];   //[ngamma]
   Int_t           ngenLep;
   Float_t         genLep_pt[2];   //[ngenLep]
   Float_t         genLep_eta[2];   //[ngenLep]
   Float_t         genLep_phi[2];   //[ngenLep]
   Float_t         genLep_mass[2];   //[ngenLep]
   Int_t           genLep_pdgId[2];   //[ngenLep]
   Float_t         genLep_charge[2];   //[ngenLep]
   Int_t           genLep_status[2];   //[ngenLep]
   Int_t           genLep_sourceId[2];   //[ngenLep]
   Float_t         evt_scale1fb;
   Float_t         evt_xsec;
   Float_t         evt_kfactor;
   Float_t         evt_filter;
   ULong64_t       evt_nEvts;
   //Int_t       evt_nEvts;
   Int_t           evt_id;
   Float_t weight_lepsf;
   Float_t weight_lepsf_UP;
   Float_t weight_lepsf_DN;
   Float_t weight_btagsf;
   Float_t weight_btagsf_UP;
   Float_t weight_btagsf_DN;
   Float_t weight_sigtrigsf;
   Float_t weight_dileptrigsf;
   Float_t weight_phottrigsf;
   Float_t weight_pu;
   Float_t weight_isr;
   Float_t weight_scales[500];
   Float_t weight_scales_UP;
   Float_t weight_scales_DN;
   Float_t weight_pdfs_UP;
   Float_t weight_pdfs_DN;
   Int_t nLHEweight;
   Int_t LHEweight_id[1000];
   Float_t LHEweight_wgt[1000];
   Float_t LHEweight_original;
   
   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_isGolden;   //!
   TBranch        *b_isSilver;   //!
   TBranch        *b_xsec;   //!
   TBranch        *b_HLT_PFMET170;   //!
   TBranch        *b_HLT_Photon90_R9Id90_HE10_IsoM;   //!
   TBranch        *b_HLT_Photon75_R9Id90_HE10_IsoM;   //!
   TBranch        *b_HLT_Photon120;   //!
   TBranch        *b_HLT_Photon75;   //!
   TBranch        *b_HLT_Photon165_HE10;   //!
   TBranch        *b_HLT_Photon120_R9Id90_HE10_IsoM;   //!
   TBranch        *b_HLT_Photon90;   //!
   TBranch        *b_HLT_PFHT350_PFMET100;   //!
   TBranch        *b_HLT_PFMETNoMu90_PFMHTNoMu90;   //!
   TBranch        *b_HLT_PFMET90_PFMHT90;   //!
   TBranch        *b_HLT_PFHT475_Prescale;   //!
   TBranch        *b_HLT_PFHT350_Prescale;   //!
   TBranch        *b_HLT_SingleMu;   //!
   TBranch        *b_HLT_SingleEl;   //!
   TBranch        *b_HLT_MuX_Ele12;   //!
   TBranch        *b_HLT_Mu8_EleX;   //!
   TBranch        *b_HLT_PFHT800;   //!
   TBranch        *b_HLT_Photon155;   //!
   TBranch        *b_HLT_PFHT900;   //!
   TBranch        *b_HLT_Photon175;   //!
   TBranch        *b_HLT_DoubleEl;   //!
   TBranch        *b_HLT_DoubleEle33;   //!
   TBranch        *b_HLT_DoubleMu;   //!
   TBranch        *b_HLT_Mu30_TkMu11;
   //TBranch        *b_Flag_EcalDeadCellTriggerPrimitiveFilter;   //!
   //TBranch        *b_Flag_trkPOG_manystripclus53X;   //!
   //TBranch        *b_Flag_ecalLaserCorrFilter;   //!
   //TBranch        *b_Flag_trkPOG_toomanystripclus53X;   //!
   //TBranch        *b_Flag_hcalLaserEventFilter;   //!
   //TBranch        *b_Flag_trkPOG_logErrorTooManyClusters;   //!
   //TBranch        *b_Flag_trkPOGFilters;   //!
   //TBranch        *b_Flag_trackingFailureFilter;   //!
   //TBranch        *b_Flag_CSCTightHaloFilter;   //!
   //TBranch        *b_Flag_HBHENoiseFilter;   //!
   //TBranch        *b_Flag_HBHENoiseIsoFilter;   //!
   //TBranch        *b_Flag_goodVertices;   //!
   //TBranch        *b_Flag_METFilters;   //!
   //TBranch        *b_Flag_eeBadScFilter;   //!
   TBranch        *b_Flag_EcalDeadCellTriggerPrimitiveFilter;   //!
   TBranch        *b_Flag_HBHENoiseIsoFilter;   //!
   TBranch        *b_Flag_goodVertices;   //!
   TBranch        *b_Flag_HBHENoiseFilter;   //!
   TBranch        *b_Flag_CSCTightHalo2015Filter;   //!
   TBranch        *b_Flag_eeBadScFilter;   //!
   TBranch        *b_Flag_globalTightHalo2016Filter;   //!
   TBranch        *b_puWeight;   //!
   TBranch        *b_nTrueInt;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_nVert;   //!
   TBranch        *b_genQScale;   //!
   TBranch        *b_met_rawPt;   //!
   TBranch        *b_met_rawPhi;   //!
   TBranch        *b_nMuons10;   //!
   TBranch        *b_nElectrons10;   //!
   TBranch        *b_nTaus20;   //!
   TBranch        *b_nGammas20;   //!
   TBranch        *b_nLepLowMT;   //! 
   TBranch        *b_nPFLep5LowMT;   //!
   TBranch        *b_nPFHad10LowMT;   //!
   TBranch        *b_met_pt;   //!
   TBranch        *b_met_eta;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_met_mass;   //!
   TBranch        *b_met_sumEt;   //!
   TBranch        *b_met_genPt;   //!
   TBranch        *b_met_genPhi;   //!
   TBranch        *b_met_genEta;   //!
   TBranch        *b_nGenPart;   //!
   TBranch        *b_GenPart_motherId;   //!
   TBranch        *b_GenPart_grandmotherId;   //!
   TBranch        *b_GenPart_sourceId;   //!
   TBranch        *b_GenPart_charge;   //!
   TBranch        *b_GenPart_status;   //!
   TBranch        *b_GenPart_pdgId;   //!
   TBranch        *b_GenPart_pt;   //!
   TBranch        *b_GenPart_eta;   //!
   TBranch        *b_GenPart_phi;   //!
   TBranch        *b_GenPart_mass;   //!
   TBranch        *b_GenPart_motherIndex;   //!
   TBranch        *b_ngenPart;   //!
   TBranch        *b_genPart_pt;   //!
   TBranch        *b_genPart_eta;   //!
   TBranch        *b_genPart_phi;   //!
   TBranch        *b_genPart_mass;   //!
   TBranch        *b_genPart_pdgId;   //!
   TBranch        *b_genPart_charge;   //!
   TBranch        *b_genPart_status;   //!
   TBranch        *b_genPart_motherId;   //!
   TBranch        *b_genPart_motherIndex;   //!
   TBranch        *b_genPart_grandmaId;   //!
   TBranch        *b_njet;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_mass;   //!
   TBranch        *b_jet_btagCSV;   //!
   TBranch        *b_jet_rawPt;   //!
   TBranch        *b_jet_mcPt;   //!
   TBranch        *b_jet_mcFlavour;   //!
   TBranch        *b_jet_mcMatchId;   //!
   TBranch        *b_jet_mcMatchFlav;   //!
   TBranch        *b_jet_area;   //!
   TBranch        *b_jet_puId;   //!
   TBranch        *b_jet_id;   //!
   TBranch        *b_jet_ptd;   //!
   TBranch        *b_jet_axis2;   //!
   TBranch        *b_jet_mult;   //!
   TBranch        *b_jet_qgl;   //!
   TBranch        *b_jet_partonId;   //!
   TBranch        *b_jet_partonMotherId;   //!
   TBranch        *b_jet_chHEF;   //!
   TBranch        *b_jet_neHEF;   //!
   TBranch        *b_jet_phEF;   //!
   TBranch        *b_jet_eEF;   //!
   TBranch        *b_jet_muEF;   //!
   TBranch        *b_jet_HFHEF;   //!
   TBranch        *b_jet_HFEMEF;   //!
   TBranch        *b_jet_chHMult;   //!
   TBranch        *b_jet_neHMult;   //!
   TBranch        *b_jet_phMult;   //!
   TBranch        *b_jet_eMult;   //!
   TBranch        *b_jet_muMult;   //!
   TBranch        *b_jet_HFHMult;   //!
   TBranch        *b_jet_HFEMMult;   //!
   TBranch        *b_nlep;   //!
   TBranch        *b_lep_pt;   //!
   TBranch        *b_lep_eta;   //!
   TBranch        *b_lep_phi;   //!
   TBranch        *b_lep_mass;   //!
   TBranch        *b_lep_pdgId;   //!
   TBranch        *b_lep_charge;   //!
   TBranch        *b_lep_dxy;   //!
   TBranch        *b_lep_dz;   //!
   TBranch        *b_lep_edxy;   //!
   TBranch        *b_lep_edz;   //!
   TBranch        *b_lep_ip3d;   //!
   TBranch        *b_lep_sip3d;   //!
   TBranch        *b_lep_tightId;   //!
   TBranch        *b_lep_convVeto;   //!
   TBranch        *b_lep_lostHits;   //!
   TBranch        *b_lep_looseIdSusy;   //!
   TBranch        *b_lep_relIso03;   //!
   TBranch        *b_lep_relIso04;   //!
   TBranch        *b_lep_miniRelIso;   //!
   TBranch        *b_lep_chargedHadRelIso03;   //!
   TBranch        *b_lep_chargedHadRelIso04;   //!
   TBranch        *b_lep_convVetoFull;   //!
   TBranch        *b_lep_eleCutIdCSA14_25ns_v1;   //!
   TBranch        *b_lep_eleCutIdCSA14_50ns_v1;   //!
   TBranch        *b_lep_eleMVAId;   //!
   TBranch        *b_lep_tightCharge;   //!
   TBranch        *b_lep_mvaId;   //!
   TBranch        *b_lep_mvaIdTrig;   //!
   TBranch        *b_lep_nStations;   //!
   TBranch        *b_lep_trkKink;   //!
   TBranch        *b_lep_caloCompatibility;   //!
   TBranch        *b_lep_globalTrackChi2;   //!
   TBranch        *b_lep_trackerLayers;   //!
   TBranch        *b_lep_pixelLayers;   //!
   TBranch        *b_lep_mvaTTH;   //!
   TBranch        *b_lep_mvaSusy;   //!
   TBranch        *b_lep_jetPtRatio;   //!
   TBranch        *b_lep_jetBTagCSV;   //!
   TBranch        *b_lep_jetDR;   //!
   TBranch        *b_lep_mcMatchId;   //!
   TBranch        *b_lep_mcMatchAny;   //!
   TBranch        *b_lep_mcMatchTau;   //!
   TBranch        *b_lep_r9;   //!
   TBranch        *b_ngenLepFromTau;   //!
   TBranch        *b_genLepFromTau_pt;   //!
   TBranch        *b_genLepFromTau_eta;   //!
   TBranch        *b_genLepFromTau_phi;   //!
   TBranch        *b_genLepFromTau_mass;   //!
   TBranch        *b_genLepFromTau_pdgId;   //!
   TBranch        *b_genLepFromTau_charge;   //!
   TBranch        *b_genLepFromTau_status;   //!
   TBranch        *b_genLepFromTau_sourceId;   //!
   TBranch        *b_ngamma;   //!
   TBranch        *b_gamma_pt;   //!
   TBranch        *b_gamma_eta;   //!
   TBranch        *b_gamma_phi;   //!
   TBranch        *b_gamma_mass;   //!
   TBranch        *b_gamma_pdgId;   //!
   TBranch        *b_gamma_idCutBased;   //!
   TBranch        *b_gamma_hOverE;   //!
   TBranch        *b_gamma_r9;   //!
   TBranch        *b_gamma_sigmaIetaIeta;   //!
   TBranch        *b_gamma_chHadIso;   //!
   TBranch        *b_gamma_chHadIso04;   //!
   TBranch        *b_gamma_chHadIsoRC;   //!
   TBranch        *b_gamma_chHadIsoRC04;   //!
   TBranch        *b_gamma_neuHadIso;   //!
   TBranch        *b_gamma_phIso;   //!
   TBranch        *b_gamma_mcMatchId;   //!
   TBranch        *b_gamma_mcPt;   //!
   TBranch        *b_gamma_genIso04;   //!
   TBranch        *b_gamma_genIso03;   //!
   TBranch        *b_gamma_drMinParton;   //!
   TBranch        *b_ngenLep;   //!
   TBranch        *b_genLep_pt;   //!
   TBranch        *b_genLep_eta;   //!
   TBranch        *b_genLep_phi;   //!
   TBranch        *b_genLep_mass;   //!
   TBranch        *b_genLep_pdgId;   //!
   TBranch        *b_genLep_charge;   //!
   TBranch        *b_genLep_status;   //!
   TBranch        *b_genLep_sourceId;   //!
   TBranch        *b_evt_scale1fb;   //!
   TBranch        *b_evt_xsec;   //!
   TBranch        *b_evt_kfactor;   //!
   TBranch        *b_evt_filter;   //!
   TBranch        *b_evt_nEvts;   //!
   TBranch        *b_evt_id;   //!
   TBranch *b_weight_lepsf;
   TBranch *b_weight_lepsf_UP;
   TBranch *b_weight_lepsf_DN;
   TBranch *b_weight_btagsf;
   TBranch *b_weight_btagsf_UP;
   TBranch *b_weight_btagsf_DN;
   TBranch *b_weight_sigtrigsf;
   TBranch *b_weight_dileptrigsf;
   TBranch *b_weight_phottrigsf;
   TBranch *b_weight_pu;
   TBranch *b_weight_isr;
   TBranch *b_weight_scales;
   TBranch *b_weight_scales_UP;
   TBranch *b_weight_scales_DN;
   TBranch *b_weight_pdfs_UP;
   TBranch *b_weight_pdfs_DN;
   TBranch *b_nLHEweight;
   TBranch *b_LHEweight_id;
   TBranch *b_LHEweight_wgt;
   TBranch *b_LHEweight_original;

   bool loadGenStuff;

   ZGTree(TTree *tree=0);
   virtual ~ZGTree();
   virtual Bool_t   passFilters     () const;

   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   //virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef zg_cxx
ZGTree::ZGTree(TTree *tree) : fChain(0) 
{

  loadGenStuff = true;

}

ZGTree::~ZGTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ZGTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ZGTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ZGTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("isGolden", &isGolden, &b_isGolden);
   fChain->SetBranchAddress("isSilver", &isSilver, &b_isSilver);
   fChain->SetBranchAddress("xsec", &xsec, &b_xsec);
   fChain->SetBranchAddress("HLT_PFMET170", &HLT_PFMET170, &b_HLT_PFMET170);
   fChain->SetBranchAddress("HLT_Photon90_R9Id90_HE10_IsoM", &HLT_Photon90_R9Id90_HE10_IsoM, &b_HLT_Photon90_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM", &HLT_Photon75_R9Id90_HE10_IsoM, &b_HLT_Photon75_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon120", &HLT_Photon120, &b_HLT_Photon120);
   fChain->SetBranchAddress("HLT_Photon75", &HLT_Photon75, &b_HLT_Photon75);
   fChain->SetBranchAddress("HLT_Photon165_HE10", &HLT_Photon165_HE10, &b_HLT_Photon165_HE10);
   fChain->SetBranchAddress("HLT_Photon120_R9Id90_HE10_IsoM", &HLT_Photon120_R9Id90_HE10_IsoM, &b_HLT_Photon120_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon90", &HLT_Photon90, &b_HLT_Photon90);
   fChain->SetBranchAddress("HLT_PFHT350_PFMET100", &HLT_PFHT350_PFMET100, &b_HLT_PFHT350_PFMET100);
   fChain->SetBranchAddress("HLT_PFMETNoMu90_PFMHTNoMu90", &HLT_PFMETNoMu90_PFMHTNoMu90, &b_HLT_PFMETNoMu90_PFMHTNoMu90);
   fChain->SetBranchAddress("HLT_PFMET90_PFMHT90", &HLT_PFMET90_PFMHT90, &b_HLT_PFMET90_PFMHT90);
   fChain->SetBranchAddress("HLT_PFHT475_Prescale", &HLT_PFHT475_Prescale, &b_HLT_PFHT475_Prescale);
   fChain->SetBranchAddress("HLT_PFHT350_Prescale", &HLT_PFHT350_Prescale, &b_HLT_PFHT350_Prescale);
   fChain->SetBranchAddress("HLT_SingleMu", &HLT_SingleMu, &b_HLT_SingleMu);
   fChain->SetBranchAddress("HLT_SingleEl", &HLT_SingleEl, &b_HLT_SingleEl);
   fChain->SetBranchAddress("HLT_Mu8_EleX", &HLT_Mu8_EleX, &b_HLT_Mu8_EleX);
   fChain->SetBranchAddress("HLT_MuX_Ele12", &HLT_MuX_Ele12, &b_HLT_MuX_Ele12);
   fChain->SetBranchAddress("HLT_PFHT800", &HLT_PFHT800, &b_HLT_PFHT800);
   fChain->SetBranchAddress("HLT_Photon155", &HLT_Photon155, &b_HLT_Photon155);
   fChain->SetBranchAddress("HLT_PFHT900", &HLT_PFHT900, &b_HLT_PFHT900);
   fChain->SetBranchAddress("HLT_Photon175", &HLT_Photon175, &b_HLT_Photon175);
   fChain->SetBranchAddress("HLT_DoubleEl", &HLT_DoubleEl, &b_HLT_DoubleEl);
   fChain->SetBranchAddress("HLT_DoubleEle33", &HLT_DoubleEle33, &b_HLT_DoubleEle33);
   fChain->SetBranchAddress("HLT_DoubleMu", &HLT_DoubleMu, &b_HLT_DoubleMu);
   fChain->SetBranchAddress("HLT_Mu30_TkMu11", &HLT_Mu30_TkMu11, &b_HLT_Mu30_TkMu11);
   //fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
   //fChain->SetBranchAddress("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X, &b_Flag_trkPOG_manystripclus53X);
   //fChain->SetBranchAddress("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter, &b_Flag_ecalLaserCorrFilter);
   //fChain->SetBranchAddress("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X, &b_Flag_trkPOG_toomanystripclus53X);
   //fChain->SetBranchAddress("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter, &b_Flag_hcalLaserEventFilter);
   //fChain->SetBranchAddress("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters, &b_Flag_trkPOG_logErrorTooManyClusters);
   //fChain->SetBranchAddress("Flag_trkPOGFilters", &Flag_trkPOGFilters, &b_Flag_trkPOGFilters);
   //fChain->SetBranchAddress("Flag_trackingFailureFilter", &Flag_trackingFailureFilter, &b_Flag_trackingFailureFilter);
   //fChain->SetBranchAddress("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter, &b_Flag_CSCTightHaloFilter);
   //fChain->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
   //fChain->SetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter, &b_Flag_HBHENoiseIsoFilter);
   //fChain->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
   //fChain->SetBranchAddress("Flag_METFilters", &Flag_METFilters, &b_Flag_METFilters);
   //fChain->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter, &b_Flag_HBHENoiseIsoFilter);
   fChain->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
   fChain->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
   fChain->SetBranchAddress("Flag_CSCTightHalo2015Filter", &Flag_CSCTightHalo2015Filter, &b_Flag_CSCTightHalo2015Filter);
   fChain->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
   fChain->SetBranchAddress("Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter, &b_Flag_globalTightHalo2016Filter);
   fChain->SetBranchAddress("puWeight", &puWeight, &b_puWeight);
   fChain->SetBranchAddress("nTrueInt", &nTrueInt, &b_nTrueInt);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("nVert", &nVert, &b_nVert);
   fChain->SetBranchAddress("genQScale", &genQScale, &b_genQScale);
   fChain->SetBranchAddress("met_rawPt", &met_rawPt, &b_met_rawPt);
   fChain->SetBranchAddress("met_rawPhi", &met_rawPhi, &b_met_rawPhi);
   fChain->SetBranchAddress("nMuons10", &nMuons10, &b_nMuons10);
   fChain->SetBranchAddress("nElectrons10", &nElectrons10, &b_nElectrons10);
   fChain->SetBranchAddress("nTaus20", &nTaus20, &b_nTaus20);
   fChain->SetBranchAddress("nGammas20", &nGammas20, &b_nGammas20);
   fChain->SetBranchAddress("nLepLowMT", &nLepLowMT, &b_nLepLowMT);
   fChain->SetBranchAddress("nPFLep5LowMT", &nPFLep5LowMT, &b_nPFLep5LowMT);
   fChain->SetBranchAddress("nPFHad10LowMT", &nPFHad10LowMT, &b_nPFHad10LowMT);
   if( loadGenStuff ) {
     fChain->SetBranchAddress("nGenPart", &nGenPart, &b_nGenPart);
     fChain->SetBranchAddress("GenPart_motherId", GenPart_motherId, &b_GenPart_motherId);
     fChain->SetBranchAddress("GenPart_grandmotherId", GenPart_grandmotherId, &b_GenPart_grandmotherId);
     fChain->SetBranchAddress("GenPart_sourceId", GenPart_sourceId, &b_GenPart_sourceId);
     fChain->SetBranchAddress("GenPart_charge", GenPart_charge, &b_GenPart_charge);
     fChain->SetBranchAddress("GenPart_status", GenPart_status, &b_GenPart_status);
     fChain->SetBranchAddress("GenPart_pdgId", GenPart_pdgId, &b_GenPart_pdgId);
     fChain->SetBranchAddress("GenPart_pt", GenPart_pt, &b_GenPart_pt);
     fChain->SetBranchAddress("GenPart_eta", GenPart_eta, &b_GenPart_eta);
     fChain->SetBranchAddress("GenPart_phi", GenPart_phi, &b_GenPart_phi);
     fChain->SetBranchAddress("GenPart_mass", GenPart_mass, &b_GenPart_mass);
     fChain->SetBranchAddress("GenPart_motherIndex", GenPart_motherIndex, &b_GenPart_motherIndex);
     fChain->SetBranchAddress("ngenPart", &ngenPart, &b_ngenPart);
     fChain->SetBranchAddress("genPart_pt", genPart_pt, &b_genPart_pt);
     fChain->SetBranchAddress("genPart_eta", genPart_eta, &b_genPart_eta);
     fChain->SetBranchAddress("genPart_phi", genPart_phi, &b_genPart_phi);
     fChain->SetBranchAddress("genPart_mass", genPart_mass, &b_genPart_mass);
     fChain->SetBranchAddress("genPart_pdgId", genPart_pdgId, &b_genPart_pdgId);
     fChain->SetBranchAddress("genPart_charge", genPart_charge, &b_genPart_charge);
     fChain->SetBranchAddress("genPart_status", genPart_status, &b_genPart_status);
     fChain->SetBranchAddress("genPart_motherId", genPart_motherId, &b_genPart_motherId);
     fChain->SetBranchAddress("genPart_motherIndex", genPart_motherIndex, &b_genPart_motherIndex);
     fChain->SetBranchAddress("genPart_grandmaId", genPart_grandmaId, &b_genPart_grandmaId);
     fChain->SetBranchAddress("ngenLepFromTau", &ngenLepFromTau, &b_ngenLepFromTau);
     fChain->SetBranchAddress("genLepFromTau_pt", genLepFromTau_pt, &b_genLepFromTau_pt);
     fChain->SetBranchAddress("genLepFromTau_eta", genLepFromTau_eta, &b_genLepFromTau_eta);
     fChain->SetBranchAddress("genLepFromTau_phi", genLepFromTau_phi, &b_genLepFromTau_phi);
     fChain->SetBranchAddress("genLepFromTau_mass", genLepFromTau_mass, &b_genLepFromTau_mass);
     fChain->SetBranchAddress("genLepFromTau_pdgId", genLepFromTau_pdgId, &b_genLepFromTau_pdgId);
     fChain->SetBranchAddress("genLepFromTau_charge", genLepFromTau_charge, &b_genLepFromTau_charge);
     fChain->SetBranchAddress("genLepFromTau_status", genLepFromTau_status, &b_genLepFromTau_status);
     fChain->SetBranchAddress("genLepFromTau_sourceId", genLepFromTau_sourceId, &b_genLepFromTau_sourceId);
     fChain->SetBranchAddress("ngenLep", &ngenLep, &b_ngenLep);
     fChain->SetBranchAddress("genLep_pt", genLep_pt, &b_genLep_pt);
     fChain->SetBranchAddress("genLep_eta", genLep_eta, &b_genLep_eta);
     fChain->SetBranchAddress("genLep_phi", genLep_phi, &b_genLep_phi);
     fChain->SetBranchAddress("genLep_mass", genLep_mass, &b_genLep_mass);
     fChain->SetBranchAddress("genLep_pdgId", genLep_pdgId, &b_genLep_pdgId);
     fChain->SetBranchAddress("genLep_charge", genLep_charge, &b_genLep_charge);
     fChain->SetBranchAddress("genLep_status", genLep_status, &b_genLep_status);
     fChain->SetBranchAddress("genLep_sourceId", genLep_sourceId, &b_genLep_sourceId);
   }
   fChain->SetBranchAddress("njet", &njet, &b_njet);
   fChain->SetBranchAddress("jet_pt", jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_eta", jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_mass", jet_mass, &b_jet_mass);
   fChain->SetBranchAddress("jet_btagCSV", jet_btagCSV, &b_jet_btagCSV);
   fChain->SetBranchAddress("jet_rawPt", jet_rawPt, &b_jet_rawPt);
   fChain->SetBranchAddress("jet_mcPt", jet_mcPt, &b_jet_mcPt);
   fChain->SetBranchAddress("jet_mcFlavour", jet_mcFlavour, &b_jet_mcFlavour);
   fChain->SetBranchAddress("jet_mcMatchId", jet_mcMatchId, &b_jet_mcMatchId);
   fChain->SetBranchAddress("jet_mcMatchFlav", jet_mcMatchFlav, &b_jet_mcMatchFlav);
   fChain->SetBranchAddress("jet_area", jet_area, &b_jet_area);
   fChain->SetBranchAddress("jet_puId", jet_puId, &b_jet_puId);
   fChain->SetBranchAddress("jet_id", jet_id, &b_jet_id);
   fChain->SetBranchAddress("jet_ptd", jet_ptd, &b_jet_ptd);
   fChain->SetBranchAddress("jet_axis2", jet_axis2, &b_jet_axis2);
   fChain->SetBranchAddress("jet_mult", jet_mult, &b_jet_mult);
   fChain->SetBranchAddress("jet_qgl", jet_qgl, &b_jet_qgl);
   fChain->SetBranchAddress("jet_partonId", jet_partonId, &b_jet_partonId);
   fChain->SetBranchAddress("jet_partonMotherId", jet_partonMotherId, &b_jet_partonMotherId);
   fChain->SetBranchAddress("jet_chHEF", jet_chHEF, &b_jet_chHEF);
   fChain->SetBranchAddress("jet_neHEF", jet_neHEF, &b_jet_neHEF);
   fChain->SetBranchAddress("jet_phEF", jet_phEF, &b_jet_phEF);
   fChain->SetBranchAddress("jet_eEF", jet_eEF, &b_jet_eEF);
   fChain->SetBranchAddress("jet_muEF", jet_muEF, &b_jet_muEF);
   fChain->SetBranchAddress("jet_HFHEF", jet_HFHEF, &b_jet_HFHEF);
   fChain->SetBranchAddress("jet_HFEMEF", jet_HFEMEF, &b_jet_HFEMEF);
   fChain->SetBranchAddress("jet_chHMult", jet_chHMult, &b_jet_chHMult);
   fChain->SetBranchAddress("jet_neHMult", jet_neHMult, &b_jet_neHMult);
   fChain->SetBranchAddress("jet_phMult", jet_phMult, &b_jet_phMult);
   fChain->SetBranchAddress("jet_eMult", jet_eMult, &b_jet_eMult);
   fChain->SetBranchAddress("jet_muMult", jet_muMult, &b_jet_muMult);
   fChain->SetBranchAddress("jet_HFHMult", jet_HFHMult, &b_jet_HFHMult);
   fChain->SetBranchAddress("jet_HFEMMult", jet_HFEMMult, &b_jet_HFEMMult);
   fChain->SetBranchAddress("nlep", &nlep, &b_nlep);
   fChain->SetBranchAddress("lep_pt", lep_pt, &b_lep_pt);
   fChain->SetBranchAddress("lep_eta", lep_eta, &b_lep_eta);
   fChain->SetBranchAddress("lep_phi", lep_phi, &b_lep_phi);
   fChain->SetBranchAddress("lep_mass", lep_mass, &b_lep_mass);
   fChain->SetBranchAddress("lep_pdgId", lep_pdgId, &b_lep_pdgId);
   fChain->SetBranchAddress("lep_charge", lep_charge, &b_lep_charge);
   fChain->SetBranchAddress("lep_dxy", lep_dxy, &b_lep_dxy);
   fChain->SetBranchAddress("lep_dz", lep_dz, &b_lep_dz);
   fChain->SetBranchAddress("lep_edxy", lep_edxy, &b_lep_edxy);
   fChain->SetBranchAddress("lep_edz", lep_edz, &b_lep_edz);
   fChain->SetBranchAddress("lep_ip3d", lep_ip3d, &b_lep_ip3d);
   fChain->SetBranchAddress("lep_sip3d", lep_sip3d, &b_lep_sip3d);
   fChain->SetBranchAddress("lep_tightId", lep_tightId, &b_lep_tightId);
   fChain->SetBranchAddress("lep_convVeto", lep_convVeto, &b_lep_convVeto);
   fChain->SetBranchAddress("lep_lostHits", lep_lostHits, &b_lep_lostHits);
   fChain->SetBranchAddress("lep_looseIdSusy", lep_looseIdSusy, &b_lep_looseIdSusy);
   fChain->SetBranchAddress("lep_relIso03", lep_relIso03, &b_lep_relIso03);
   fChain->SetBranchAddress("lep_relIso04", lep_relIso04, &b_lep_relIso04);
   fChain->SetBranchAddress("lep_miniRelIso", lep_miniRelIso, &b_lep_miniRelIso);
   fChain->SetBranchAddress("lep_chargedHadRelIso03", lep_chargedHadRelIso03, &b_lep_chargedHadRelIso03);
   fChain->SetBranchAddress("lep_chargedHadRelIso04", lep_chargedHadRelIso04, &b_lep_chargedHadRelIso04);
   fChain->SetBranchAddress("lep_convVetoFull", lep_convVetoFull, &b_lep_convVetoFull);
   fChain->SetBranchAddress("lep_eleCutIdCSA14_25ns_v1", lep_eleCutIdCSA14_25ns_v1, &b_lep_eleCutIdCSA14_25ns_v1);
   fChain->SetBranchAddress("lep_eleCutIdCSA14_50ns_v1", lep_eleCutIdCSA14_50ns_v1, &b_lep_eleCutIdCSA14_50ns_v1);
   fChain->SetBranchAddress("lep_eleMVAId", lep_eleMVAId, &b_lep_eleMVAId);
   fChain->SetBranchAddress("lep_tightCharge", lep_tightCharge, &b_lep_tightCharge);
   fChain->SetBranchAddress("lep_mvaId", lep_mvaId, &b_lep_mvaId);
   fChain->SetBranchAddress("lep_mvaIdTrig", lep_mvaIdTrig, &b_lep_mvaIdTrig);
   fChain->SetBranchAddress("lep_nStations", lep_nStations, &b_lep_nStations);
   fChain->SetBranchAddress("lep_trkKink", lep_trkKink, &b_lep_trkKink);
   fChain->SetBranchAddress("lep_caloCompatibility", lep_caloCompatibility, &b_lep_caloCompatibility);
   fChain->SetBranchAddress("lep_globalTrackChi2", lep_globalTrackChi2, &b_lep_globalTrackChi2);
   fChain->SetBranchAddress("lep_trackerLayers", lep_trackerLayers, &b_lep_trackerLayers);
   fChain->SetBranchAddress("lep_pixelLayers", lep_pixelLayers, &b_lep_pixelLayers);
   fChain->SetBranchAddress("lep_mvaTTH", lep_mvaTTH, &b_lep_mvaTTH);
   fChain->SetBranchAddress("lep_mvaSusy", lep_mvaSusy, &b_lep_mvaSusy);
   fChain->SetBranchAddress("lep_jetPtRatio", lep_jetPtRatio, &b_lep_jetPtRatio);
   fChain->SetBranchAddress("lep_jetBTagCSV", lep_jetBTagCSV, &b_lep_jetBTagCSV);
   fChain->SetBranchAddress("lep_jetDR", lep_jetDR, &b_lep_jetDR);
   fChain->SetBranchAddress("lep_mcMatchId", lep_mcMatchId, &b_lep_mcMatchId);
   fChain->SetBranchAddress("lep_mcMatchAny", lep_mcMatchAny, &b_lep_mcMatchAny);
   fChain->SetBranchAddress("lep_mcMatchTau", lep_mcMatchTau, &b_lep_mcMatchTau);
   fChain->SetBranchAddress("lep_r9", lep_r9, &b_lep_r9);
   fChain->SetBranchAddress("ngamma", &ngamma, &b_ngamma);
   fChain->SetBranchAddress("gamma_pt", gamma_pt, &b_gamma_pt);
   fChain->SetBranchAddress("gamma_eta", gamma_eta, &b_gamma_eta);
   fChain->SetBranchAddress("gamma_phi", gamma_phi, &b_gamma_phi);
   fChain->SetBranchAddress("gamma_mass", gamma_mass, &b_gamma_mass);
   fChain->SetBranchAddress("gamma_pdgId", gamma_pdgId, &b_gamma_pdgId);
   fChain->SetBranchAddress("gamma_idCutBased", gamma_idCutBased, &b_gamma_idCutBased);
   fChain->SetBranchAddress("gamma_hOverE", gamma_hOverE, &b_gamma_hOverE);
   fChain->SetBranchAddress("gamma_r9", gamma_r9, &b_gamma_r9);
   fChain->SetBranchAddress("gamma_sigmaIetaIeta", gamma_sigmaIetaIeta, &b_gamma_sigmaIetaIeta);
   fChain->SetBranchAddress("gamma_chHadIso", gamma_chHadIso, &b_gamma_chHadIso);
   fChain->SetBranchAddress("gamma_chHadIso04", gamma_chHadIso04, &b_gamma_chHadIso04);
   fChain->SetBranchAddress("gamma_chHadIsoRC"  , gamma_chHadIsoRC,   &b_gamma_chHadIsoRC);
   fChain->SetBranchAddress("gamma_chHadIsoRC04", gamma_chHadIsoRC04, &b_gamma_chHadIsoRC04);
   fChain->SetBranchAddress("gamma_neuHadIso", gamma_neuHadIso, &b_gamma_neuHadIso);
   fChain->SetBranchAddress("gamma_phIso", gamma_phIso, &b_gamma_phIso);
   fChain->SetBranchAddress("gamma_mcMatchId", gamma_mcMatchId, &b_gamma_mcMatchId);
   fChain->SetBranchAddress("gamma_mcPt", gamma_mcPt, &b_gamma_mcPt);
   fChain->SetBranchAddress("gamma_genIso04", gamma_genIso04, &b_gamma_genIso04);
   fChain->SetBranchAddress("gamma_genIso03", gamma_genIso03, &b_gamma_genIso03);
   fChain->SetBranchAddress("gamma_drMinParton", gamma_drMinParton, &b_gamma_drMinParton);
   fChain->SetBranchAddress("evt_scale1fb", &evt_scale1fb, &b_evt_scale1fb);
   fChain->SetBranchAddress("evt_xsec", &evt_xsec, &b_evt_xsec);
   fChain->SetBranchAddress("evt_kfactor", &evt_kfactor, &b_evt_kfactor);
   fChain->SetBranchAddress("evt_filter", &evt_filter, &b_evt_filter);
   fChain->SetBranchAddress("evt_nEvts", &evt_nEvts, &b_evt_nEvts);
   fChain->SetBranchAddress("evt_id", &evt_id, &b_evt_id);
   fChain->SetBranchAddress("weight_lepsf", &weight_lepsf, &b_weight_lepsf);
   fChain->SetBranchAddress("weight_lepsf_UP", &weight_lepsf_UP, &b_weight_lepsf_UP);
   fChain->SetBranchAddress("weight_lepsf_DN", &weight_lepsf_DN, &b_weight_lepsf_DN);
   fChain->SetBranchAddress("weight_btagsf", &weight_btagsf, &b_weight_btagsf);
   fChain->SetBranchAddress("weight_btagsf_UP", &weight_btagsf_UP, &b_weight_btagsf_UP);
   fChain->SetBranchAddress("weight_btagsf_DN", &weight_btagsf_DN, &b_weight_btagsf_DN);
   fChain->SetBranchAddress("weight_sigtrigsf", &weight_sigtrigsf, &b_weight_sigtrigsf);
   fChain->SetBranchAddress("weight_dileptrigsf", &weight_dileptrigsf, &b_weight_dileptrigsf);
   fChain->SetBranchAddress("weight_phottrigsf", &weight_phottrigsf, &b_weight_phottrigsf);
   fChain->SetBranchAddress("weight_pu", &weight_pu, &b_weight_pu);
   fChain->SetBranchAddress("weight_isr", &weight_isr, &b_weight_isr);
   fChain->SetBranchAddress("weight_scales", weight_scales, &b_weight_scales);
   fChain->SetBranchAddress("weight_scales_UP", &weight_scales_UP, &b_weight_scales_UP);
   fChain->SetBranchAddress("weight_scales_DN", &weight_scales_DN, &b_weight_scales_DN);
   fChain->SetBranchAddress("weight_pdfs_UP", &weight_pdfs_UP, &b_weight_pdfs_UP);
   fChain->SetBranchAddress("weight_pdfs_DN", &weight_pdfs_UP, &b_weight_pdfs_UP);
   fChain->SetBranchAddress("nLHEweight", &nLHEweight, &b_nLHEweight);
   fChain->SetBranchAddress("LHEweight_id", LHEweight_id, &b_LHEweight_id);
   fChain->SetBranchAddress("LHEweight_wgt", LHEweight_wgt, &b_LHEweight_wgt);
   fChain->SetBranchAddress("LHEweight_original", &LHEweight_original, &b_LHEweight_original);

   Notify();
}

Bool_t ZGTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ZGTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}



Bool_t ZGTree::passFilters() const {
  bool pass = Flag_HBHENoiseFilter 
           && Flag_HBHENoiseIsoFilter 
           && Flag_CSCTightHalo2015Filter 
           && Flag_EcalDeadCellTriggerPrimitiveFilter 
           && Flag_goodVertices 
           && Flag_eeBadScFilter;
  return pass;
  //return nVert>0 && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_eeBadScFilter; // Beam halo from txt file
}






Int_t ZGTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef mt2_cxx
