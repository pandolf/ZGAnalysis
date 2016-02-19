#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TLegend.h"

#include "../interface/ZGDrawTools.h"




TF1* scaleFunction( TF1* f1_acc, TF1* f1_eff );

TF1* f1_acc_0p014, *f1_acc_1p4, *f1_acc_5p6;
TF1* f1_eff, *f1_eff_ee, *f1_eff_mm;


double f1_prod_0p014(double *x, double *par) 
{
    return TMath::Abs(f1_acc_0p014->EvalPar(x,par) * f1_eff->EvalPar(x,par));
}

double f1_prod_0p014_ee(double *x, double *par) 
{
    return TMath::Abs(f1_acc_0p014->EvalPar(x,par) * f1_eff_ee->EvalPar(x,par));
}

double f1_prod_0p014_mm(double *x, double *par) 
{
    return TMath::Abs(f1_acc_0p014->EvalPar(x,par) * f1_eff_mm->EvalPar(x,par));
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

  std::string file_acc_name = "genAcceptance/genAcceptance.root";
  std::cout << "-> Reading gen acceptance from: " << file_acc_name << std::endl;
  TFile* file_acc = TFile::Open(file_acc_name.c_str());

  f1_acc_0p014 = (TF1*)file_acc->Get("f1_gr_0p014");
  f1_acc_1p4   = (TF1*)file_acc->Get("f1_gr_1p4");
  f1_acc_5p6   = (TF1*)file_acc->Get("f1_gr_5p6");

  TFile* file_eff = TFile::Open("genEfficiencyFunctions.root");
  f1_eff       = (TF1*)file_eff->Get("line_all");
  f1_eff_ee    = (TF1*)file_eff->Get("line_ee");
  f1_eff_mm    = (TF1*)file_eff->Get("line_mm");

  //f1_eff = new TF1( "f1_eff", "[0]", xMin, xMax );
  //f1_eff->SetParameter(0, 0.6);

  TF1* f1_aXe_0p014    = scaleFunction( f1_acc_0p014, f1_eff    );
  TF1* f1_aXe_0p014_ee = scaleFunction( f1_acc_0p014, f1_eff_ee );
  TF1* f1_aXe_0p014_mm = scaleFunction( f1_acc_0p014, f1_eff_mm );

  //TF1* f1_aXe_0p014    = new TF1( "aXe_0p014"   , f1_prod_0p014   , xMin, xMax, 0 );
  //TF1* f1_aXe_0p014_ee = new TF1( "aXe_0p014_ee", f1_prod_0p014_ee, xMin, xMax, 0 );
  //TF1* f1_aXe_0p014_mm = new TF1( "aXe_0p014_mm", f1_prod_0p014_mm, xMin, xMax, 0 );

  TF1* f1_aXe_1p4 = scaleFunction( f1_acc_1p4, f1_eff );
  TF1* f1_aXe_5p6 = scaleFunction( f1_acc_5p6, f1_eff );


  //f1_aXe_0p014   ->SetName("aXe_0p014"   );
  //f1_aXe_0p014_ee->SetName("aXe_0p014_ee");
  //f1_aXe_0p014_mm->SetName("aXe_0p014_mm");

  //f1_aXe_1p4  ->SetName("aXe_1p4"  );
  //f1_aXe_5p6  ->SetName("aXe_5p6"  );
  

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();
  
  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, 6500, 10, 0., 1.0001 );
  h2_axes->SetXTitle( "Mass [GeV]" );
  h2_axes->SetYTitle( "Acceptance #times Efficiency" );
  h2_axes->Draw("");

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

  c1->Clear();
  c1->SetLogx(false);

  h2_axes->Draw();
  
  f1_aXe_0p014_ee->SetLineColor(46);
  f1_aXe_0p014_mm->SetLineColor(38);

  f1_aXe_0p014_ee->SetLineWidth(3);
  f1_aXe_0p014_mm->SetLineWidth(3);

  f1_aXe_0p014_ee->Draw("same");
  f1_aXe_0p014_mm->Draw("same");

  TLegend* legend2 = new TLegend( 0.2, 0.7, 0.5, 0.9 );
  legend2->SetFillColor( 0 );
  legend2->SetTextFont( 42 );
  legend2->SetHeader( "W = 0.014\%" );
  legend2->SetTextSize( 0.038 );
  legend2->AddEntry( f1_aXe_0p014_ee, "e^{+}e^{-}#gamma", "L" );
  legend2->AddEntry( f1_aXe_0p014_mm, "#mu^{+}#mu^{-}#gamma", "L" );
  legend2->Draw("same");

  ZGDrawTools::addLabels( c1, -1., "CMS Simulation" );

  gPad->RedrawAxis();

  c1->SaveAs("genAcceptanceTimesEfficiency_ee_mm.eps");
  c1->SaveAs("genAcceptanceTimesEfficiency_ee_mm.pdf");

  c1->SetLogx();
  h2_axes->GetXaxis()->SetNoExponent();
  h2_axes->GetXaxis()->SetMoreLogLabels();
  c1->SaveAs( "genAcceptanceTimesEfficiency_ee_mm_logx.eps" );
  c1->SaveAs( "genAcceptanceTimesEfficiency_ee_mm_logx.pdf" );


  TFile* file = TFile::Open("genAcceptanceTimesEfficiency.root", "recreate");
  file->cd();

  f1_aXe_0p014   ->Write();
  f1_aXe_0p014_ee->Write();
  f1_aXe_0p014_mm->Write();
  f1_aXe_1p4  ->Write();
  f1_aXe_5p6  ->Write();

  file->Close();
  
  std::cout << "-> Saved function in file: " << file->GetName() << std::endl;

  return 0;

}
  

TF1* scaleFunction( TF1* f1_acc, TF1* f1_eff ) {  // supports only scaling by const for now

  Double_t xmin, xmax;
  f1_acc->GetRange( xmin, xmax );
  TF1* f1_new = new TF1( Form("%s_times_%s", f1_acc->GetName(), f1_eff->GetName()), f1_acc->GetExpFormula(), xmin, xmax );
  for( int ipar=0; ipar<f1_acc->GetNpar(); ++ipar ) {
    f1_new->SetParameter( ipar, f1_acc->GetParameter(ipar)*f1_eff->GetParameter(0) );
  }

  return f1_new;
  
}
