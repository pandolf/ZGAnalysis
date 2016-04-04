#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLegend.h"


#include "../interface/ZGDrawTools.h"



int main() {

  ZGDrawTools::setStyle();

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  
  TH2D* h2_axes = new TH2D("axes", "", 10, 200., 2000., 10, 0., 0.0075 );
  h2_axes->SetXTitle("M(Z#gamma) [GeV]");
  h2_axes->SetYTitle("Bias [ev/GeV]");
  h2_axes->Draw();

  TF1* f1 = new TF1("Bias Function", "(x>600.)*(50000.*x^(-2.5) + 0.0001)", 200., 2000.);
  f1->SetLineColor(kRed);
  f1->Draw("same");

  ZGDrawTools::addLabels(c1, -1., "CMS Simulation" );
  gPad->RedrawAxis();

  c1->SaveAs("biasFunction.eps");
  c1->SaveAs("biasFunction.pdf");

  return 0;

}
