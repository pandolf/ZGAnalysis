#include "../interface/ZGCommonTools.h"



float ZGCommonTools::getPUweight( int nVert, TH1D* h1_data, TH1D* h1_mc ) {


  float w_data = getSinglePUweight( nVert, h1_data );
  float w_mc   = getSinglePUweight( nVert, h1_mc );

  float returnWeight = (w_mc>0.) ? w_data/w_mc : 0.;

  return returnWeight;

}


float ZGCommonTools::getSinglePUweight( int nVert, TH1D* h1 ) {

  int bin = h1->FindBin(nVert);
  float weight = h1->GetBinContent(bin)/h1->Integral("width");

  return weight;

}
