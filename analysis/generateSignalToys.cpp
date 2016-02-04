#include <iostream>

#include "RooRealVar.h"
#include "RooCBShape.h"
#include "RooDataSet.h"

#include "TFile.h"
#include "TH1D.h"



////Crystal ball function parameters are 0:alpha,1:n,2:mean,3:sigma,4:normalization;
//
//Double_t CrystalBall(Double_t *x,Double_t *par) {
//
//  Double_t t = (x[0]-par[2])/par[3];
//  if (par[0] < 0) t = -t;
//
//  Double_t absAlpha = fabs((Double_t)par[0]);
//
//  if (t >= -absAlpha) {
//    return par[4]*exp(-0.5*t*t);
//  }
//  else {
//    Double_t a =  TMath::Power(par[1]/absAlpha,par[1])*exp(-0.5*absAlpha*absAlpha);
//    Double_t b= par[1]/absAlpha - absAlpha; 
//
//    return par[4]*(a/TMath::Power(b - t, par[1]));
//  }
//}




int main() {

  std::vector<float> masses;
  masses.push_back( 300. );
  masses.push_back( 400. );
  masses.push_back( 500. );
  masses.push_back( 750. );
  masses.push_back( 1000. );
  masses.push_back( 1250. );
  masses.push_back( 1500. );
  masses.push_back( 2000. );
  masses.push_back( 2500. );
  masses.push_back( 3000. );
  masses.push_back( 5000. );


  for( unsigned i=0; i<masses.size(); ++i ) {

    std::cout << "-> Starting mass: " << masses[i] << std::endl;

    TFile* file  = TFile::Open( Form("toySignal_m%.0f.root", masses[i]), "recreate" );
    file->cd();

    float thisMass = masses[i];

    RooRealVar x("mass", "mass", thisMass, 0.2*thisMass, 2.*thisMass );

    // Crystal-Ball
    RooRealVar mean( "mean", "mean", thisMass );
    RooRealVar sigma( "sigma", "sigma", 0.015*thisMass );
    RooRealVar alpha( "alpha", "alpha", 1.2 );
    RooRealVar n( "n", "n", 0.81 );
    RooCBShape cb( "cb", "cb", x, mean, sigma, alpha, n );

    // Generate a toy dataset from the interpreted p.d.f
    RooDataSet* data = cb.generate(x,10000) ;
    data->SetName("data");
    data->Write();

    file->Close();
  
  }

  return 0;

}
