#ifndef ZGCommonTools_h
#define ZGCommonTools_h


#include "TFile.h"
#include "TH1D.h"




class ZGCommonTools {

 public:

  static float getPUweight( int nVert, TH1D* h1_data, TH1D* h1_mc );
  static float getSinglePUweight( int nVert, TH1D* h1 );

 private:

};


#endif
