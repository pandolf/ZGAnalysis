#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TLegend.h"

#include "../interface/ZGDrawTools.h"



TF1* f1_acc_0p014, *f1_acc_1p4, *f1_acc_5p6;
TF1* f1_eff;


double f1_prod_0p014(double *x, double *par) 
{
    return TMath::Abs(f1_acc_0p014->EvalPar(x,par) * f1_eff->EvalPar(x,par));
}


double f1_prod_1p4(double *x, double *par) 
{
    return TMath::Abs(f1_acc_1p4->EvalPar(x,par) * f1_eff->EvalPar(x,par));
}


double f1_prod_5p6(double *x, double *par) 
{
    return TMath::Abs(f1_acc_5p6->EvalPar(x,par) * f1_eff->EvalPar(x,par));
}



float xMin = 200.;
float xMax = 7100.;



int main() {

  ZGDrawTools::setStyle();

  TFile* file_acc = TFile::Open("genAcceptance/genAcceptance.root");

  f1_acc_0p014 = (TF1*)file_acc->Get("f1_gr_0p014");
  f1_acc_1p4   = (TF1*)file_acc->Get("f1_gr_1p4");
  f1_acc_5p6   = (TF1*)file_acc->Get("f1_gr_5p6");

  //TFile* file_eff = TFile::Open("genEfficiency.root");
  //f1_eff       = (TF1*)file_eff->Get("line_all");
  f1_eff = new TF1( "f1_eff", "[0]", xMin, xMax );
  f1_eff->SetParameter(0, 0.6);

  TF1* f1_aXe_0p014 = new TF1( "aXe_0p014", f1_prod_0p014, xMin, xMax, 0 );
  TF1* f1_aXe_1p4 = new TF1( "aXe_1p4", f1_prod_1p4, xMin, xMax, 0 );
  TF1* f1_aXe_5p6 = new TF1( "aXe_5p6", f1_prod_5p6, xMin, xMax, 0 );

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();
  
  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., 1.0001 );
  h2_axes->SetXTitle( "Mass [GeV]" );
  h2_axes->SetYTitle( "Acceptance #times Efficiency" );
  h2_axes->Draw("");

  f1_aXe_0p014->SetName("aXe_0p014");
  f1_aXe_1p4  ->SetName("aXe_1p4"  );
  f1_aXe_5p6  ->SetName("aXe_5p6"  );
  
  f1_aXe_0p014->SetLineColor(46);
  f1_aXe_1p4  ->SetLineColor(42);
  f1_aXe_5p6  ->SetLineColor(38);
  
  f1_aXe_0p014->SetLineWidth(3);
  f1_aXe_1p4  ->SetLineWidth(3);
  f1_aXe_5p6  ->SetLineWidth(3);
  
  f1_aXe_0p014->Draw("same");
  f1_aXe_1p4  ->Draw("same");
  f1_aXe_5p6  ->Draw("same");
  

  TLegend* legend = new TLegend( 0.2, 0.7, 0.5, 0.9 );
  legend->SetFillColor( 0 );
  legend->SetTextFont( 42 );
  legend->SetTextSize( 0.038 );
  legend->AddEntry( f1_aXe_0p014, "W = 0.014\%", "L" );
  legend->AddEntry( f1_aXe_1p4  , "W = 1.4\%"  , "L" );
  legend->AddEntry( f1_aXe_5p6  , "W = 5.6\%"  , "L" );
  legend->Draw("same");

  ZGDrawTools::addLabels( c1, -1., "CMS Simulation" );

  gPad->RedrawAxis();

  c1->SaveAs("genAcceptanceTimesEfficiency.eps");
  c1->SaveAs("genAcceptanceTimesEfficiency.pdf");

  c1->SetLogx();
  h2_axes->GetXaxis()->SetNoExponent();
  h2_axes->GetXaxis()->SetMoreLogLabels();
  c1->SaveAs( "genAcceptanceTimesEfficiency_logx.eps" );
  c1->SaveAs( "genAcceptanceTimesEfficiency_logx.pdf" );

  TFile* file = TFile::Open("genAcceptanceTimesEfficiency.root", "recreate");
  file->cd();

  f1_aXe_0p014->Write();
  f1_aXe_1p4  ->Write();
  f1_aXe_5p6  ->Write();

  file->Close();
  
  std::cout << "-> Saved function in file: " << file->GetName() << std::endl;

  return 0;

}
  
