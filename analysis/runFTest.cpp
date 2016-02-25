#include <iostream>
#include <fstream>

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooGenericPdf.h"
#include "RooFitResult.h"

#include "../interface/ZGConfig.h"
#include "../interface/ZGDrawTools.h"
#include "../interface/ZGCommonTools.h"


using namespace RooFit;


int doFTest( RooRealVar* x, RooDataSet* data, const std::string& funcFamily );
float fitWithModel( RooRealVar* x, RooDataSet* data, const std::string& funcFamily, unsigned iOrder );
void getPol( const std::string& name, unsigned iOrder, std::string& polFormula, RooArgSet& polargset );


int main( int argc, char* argv[] ) {


  if( argc<2 ) {
    std::cout << "USAGE: ./fitSignalShapes [configFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  std::string configFileName(argv[1]);
  ZGConfig cfg(configFileName);


  std::string eventYieldDir = cfg.getEventYieldDir();
  TFile* file = TFile::Open( Form("%s/trees.root", eventYieldDir.c_str()) );
  TTree* tree = (TTree*)file->Get("zg");

  RooRealVar* x = new RooRealVar("boss_mass", "boss_mass", 200., 200., 2200. );
  //RooRealVar* w = new RooRealVar("weight", "weight", 1., -1000., 1000. );
  RooDataSet* data = new RooDataSet( "data", "data", RooArgSet(*x), RooFit::Import(*tree) );

  std::vector<std::string> families;
  families.push_back("pow");
  families.push_back("expow");
  //families.push_back("pow");
  //families.push_back("pow");
 
  ofstream ofs("ftest.txt");
  ofs << "family     order" << std::endl;
  for( unsigned i=0; i<families.size(); ++i )
    ofs << families[i] << ": " << doFTest( x, data, families[i] ) << std::endl;

  ofs.close();

  return 0;

}



int doFTest( RooRealVar* x, RooDataSet* data, const std::string& funcFamily ) {

  ofstream ofs( Form("ftest_%s.log", funcFamily.c_str()) );

  ofs << "-> Starting F-test for family: " << funcFamily << std::endl;

  float prevNll = 0.;
  int foundOrder=-1;
  unsigned maxOrder = 5;

  for( unsigned iOrder=1; iOrder<=maxOrder && foundOrder<0; ++iOrder ) {

    float thisNll = fitWithModel( x, data, funcFamily, iOrder );

    float chi2 = 2.*( prevNll - thisNll );
    if( chi2<0. && iOrder>1 ) chi2=0.;
    float prob = TMath::Prob(chi2, 1);

    ofs << "  order: " << iOrder << " nll: " << thisNll << " chi2: " << chi2 << " Prob(chi2): " << prob << std::endl;
    if( prob>0.05 ) {
      foundOrder = iOrder;
    }

    prevNll = thisNll;

  }

  float lastGoodOrder = foundOrder-1;

  ofs << "  -> Order: " << lastGoodOrder << std::endl;
  ofs.close();

  return lastGoodOrder;

}



float fitWithModel( RooRealVar* x, RooDataSet* data, const std::string& funcFamily, unsigned iOrder ) {

  std::string polFormula;
  RooArgSet polargset;
  getPol( funcFamily, iOrder, polFormula, polargset);
  float thisNll=0.;

  RooGenericPdf* model;
  RooArgSet argset;
  argset.add( *x );
  argset.add( polargset );

  if( funcFamily=="expow" ) {

    RooRealVar* alp = new RooRealVar(Form("expow%d_alpha" , iOrder), "expow_alpha" , -4., -10., 0. );
    argset.add( *alp );
    std::string formula(Form("TMath::Max(1e-50,exp(%s)*pow(@0,@%d))", polFormula.c_str(), argset.getSize()-1));
    //std::string formula(Form("TMath::Max(1e-50,pow(@0,@%d)*exp(-(%s)))", argset.getSize()-1, polFormula.c_str()));

    std::cout << formula << std::endl;
    model = new RooGenericPdf( Form("expow_%d", iOrder), "expow", formula.c_str(), argset );

  }
  
  else if( funcFamily=="pow" ) {

    RooRealVar* alp = new RooRealVar(Form("pow%d_alpha" , iOrder), "pow_alpha" , -4., -20, 0. );

    argset.add( *alp );
    std::cout << Form("TMath::Max(1e-50,pow(%s,@%d))",polFormula.c_str(),argset.getSize()-1) << std::endl;
    model = new RooGenericPdf( Form("pow_%d", iOrder), "pow", Form("TMath::Max(1e-50,pow(%s,@%d))",polFormula.c_str(),argset.getSize()-1), argset );

  }

  RooFitResult *fitRes = model->fitTo(*data,Strategy(2),Warnings(false),Minimizer("Minuit2"),Save(true));
  fitRes->Print();
  //exit(1);
  thisNll = fitRes->minNll();

  RooPlot* frame = x->frame();
  data->plotOn(frame);
  model->plotOn(frame);
  frame->GetYaxis()->SetRangeUser( 0.01, 10000 );

  TCanvas* c1 = new TCanvas("c1", "c1", 600, 600 );
  c1->cd();
  c1->SetLogy();
  frame->Draw();

  c1->SaveAs(Form("ftest_%s_o%d.eps", funcFamily.c_str(),iOrder));

  delete c1;
  
  return thisNll;

}


void getPol( const std::string& name, unsigned iOrder, std::string& polFormula, RooArgSet& polargset ) {

  polFormula = "";

  bool coeffOnFirst = true;
  if( name=="pow" ) coeffOnFirst = false;

  for( unsigned i=1; i<=iOrder; ++i ) {

    int index = i;
    if( !coeffOnFirst ) 
      index = i-1;

    if( i!=1 ) 
      polFormula += "+";

    if( i==1 && !coeffOnFirst ) {

      polFormula = "@0";

    } else {

      std::string thisPar(Form("@%d",index));
      std::string thisTerm(thisPar);
      for( unsigned j=0; j<i; ++j )
      thisTerm += "*@0";

      polFormula += thisTerm;

      RooRealVar* par = new RooRealVar( Form("%s_o%d_p%d", name.c_str(), iOrder, index), Form("p%d", index), 0. );
      par->setConstant(false);
      //RooRealVar* par = new RooRealVar( Form("%s_o%d_p%d", name.c_str(), iOrder, index), Form("p%d", index), 0.1, 0., 20.);
      polargset.add( *par );

    }

    //if( i==1 ) { // first term has no coeff
    //  polFormula = "@0";
    //} else {
    //  polFormula += "+";
    //  std::string thisPar(Form("@%d",i-1));
    //  std::string thisTerm(thisPar);
    //  for( unsigned j=0; j<i; ++j )
    //    thisTerm += "*@0";
    //  polFormula += thisTerm;

    //  RooRealVar* par = new RooRealVar( Form("%s_o%d_p%d", name.c_str(), iOrder, i-1), Form("p%d", i-1), 0., 0., 100.);
    //  polargset.add( *par );

    //}

  } // for iorder

  std::cout << polFormula << std::endl;
  polargset.Print();

}
