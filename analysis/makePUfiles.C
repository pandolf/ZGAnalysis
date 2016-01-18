{

TChain* tree = new TChain("mt2")
tree->Add("/scratch/mmasciov/dataFromSnT_08Jan_filter/data_Run2015*_DoubleEG*.root")
TH1D* nVert = new TH1D("nVert", "", 100, 0., 100.)
tree->Project("nVert", "nVert", "HLT_DoubleEl && nlep==2 && lep_pt[1]>20.")
TFile* file= TFile::Open("puData.root", "recreate")
file->cd()
nVert->Write()
file->Close()

TFile *_file0 = TFile::Open("ZGTo2LG_post_skim.root")
TH1D* nVert = new TH1D("nVert", "", 100, 0., 100.)
nVert->Sumw2()
TTree* mt2 = (TTree*)_file0->Get("mt2");
mt2->Project("nVert", "nVert", "evt_scale1fb")
TFile* file= TFile::Open("puMC.root", "recreate")
file->cd()
nVert->Write()
file->Close()

}
