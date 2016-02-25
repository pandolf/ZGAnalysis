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


int doFTest( const std::string& outdir, RooRealVar* x, RooDataSet* data, const std::string& funcFamily, bool sumFunctions=false );
float fitWithModel( const std::string& outdir, RooRealVar* x, RooDataSet* data, const std::string& funcFamily, unsigned iOrder );
float fitWithModelSum( const std::string& outdir, RooRealVar* x, RooDataSet* data, const std::string& funcFamily, unsigned iOrder );
void getPol( const std::string& name, unsigned iOrder, std::string& polFormula, RooArgSet& polargset, float init=0., float xmin=-9999., float xmax=-9999. );
void getSumFunctionsFormula( const std::string& funcFamily, unsigned iOrder, std::string& formula, RooArgSet& argset, float init=0. );


int main( int argc, char* argv[] ) {


  if( argc<2 ) {
    std::cout << "USAGE: ./fitSignalShapes [configFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  std::string configFileName(argv[1]);
  ZGConfig cfg(configFileName);

  std::string outdir("ftest");
  system( Form("mkdir -p %s", outdir.c_str()));

  std::string eventYieldDir = cfg.getEventYieldDir();
  TFile* file = TFile::Open( Form("%s/trees.root", eventYieldDir.c_str()) );
  TTree* tree = (TTree*)file->Get("zg");

  RooRealVar* x = new RooRealVar("boss_mass", "boss_mass", 200., 200., 2200. );
  //RooRealVar* w = new RooRealVar("weight", "weight", 1., -1000., 1000. );
  RooDataSet* data = new RooDataSet( "data", "data", RooArgSet(*x), RooFit::Import(*tree) );

  std::vector<std::string> families;
  families.push_back("pow");
  families.push_back("expow");
  families.push_back("invpow");
  families.push_back("moddijet");
  families.push_back("invpowlin");
 
  std::string filename(Form("%s/ftest.txt", outdir.c_str()));
  ofstream ofs(filename.c_str());
  ofs << "family     order" << std::endl;
  for( unsigned i=0; i<families.size(); ++i )
    ofs << families[i] << ": " << doFTest( outdir, x, data, families[i] ) << std::endl;

  ofs << "sumpow: " << doFTest( outdir, x, data, "pow", true ) << std::endl;;
  ofs << "sumexp: " << doFTest( outdir, x, data, "exp", true ) << std::endl;;

  ofs.close();

  std::cout << "-> Find your results here: " << filename << std::endl;

  return 0;

}



int doFTest( const std::string& outdir, RooRealVar* x, RooDataSet* data, const std::string& funcFamily, bool sumFunctions ) {

  std::string suffix = (sumFunctions) ? "sum" : "";
  ofstream ofs( Form("%s/ftest_%s%s.log", outdir.c_str(), suffix.c_str(), funcFamily.c_str()) );

  if( sumFunctions )
    ofs << "-> Starting F-test for sum of: " << funcFamily << std::endl;
  else
    ofs << "-> Starting F-test for family: " << funcFamily << std::endl;

  float prevNll = 0.;
  int foundOrder=-1;
  unsigned maxOrder = 5;

  for( unsigned iOrder=1; iOrder<=maxOrder && foundOrder<0; ++iOrder ) {

    float thisNll = (sumFunctions) ? fitWithModelSum( outdir, x, data, funcFamily, iOrder ) : fitWithModel( outdir, x, data, funcFamily, iOrder );

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





float fitWithModel( const std::string& outdir, RooRealVar* x, RooDataSet* data, const std::string& funcFamily, unsigned iOrder ) {

  float init = 0.;
  float xmin=-99999.;
  float xmax=-99999.;
  if( funcFamily=="invpow" )
    init = 0.001;
  if( funcFamily=="moddijet" )
    init = 0.00003;

  std::string polFormula;
  RooArgSet polargset;
  getPol( funcFamily, iOrder, polFormula, polargset, init, xmin, xmax);
  float thisNll=0.;

  RooArgSet argset;
  argset.add( *x );
  argset.add( polargset );

  std::string formula;

  if( funcFamily=="expow" ) {

    RooRealVar* alp = new RooRealVar(Form("expow%d_alpha" , iOrder), "expow_alpha" , -4., -10., 0. );
    argset.add( *alp );
    formula = std::string(Form("TMath::Max(1e-50,exp(%s)*pow(@0,@%d))", polFormula.c_str(), argset.getSize()-1));

  }
  
  else if( funcFamily=="pow" ) {

    RooRealVar* alp = new RooRealVar(Form("pow%d_alpha" , iOrder), "pow_alpha" , -4., -20, 0. );
    argset.add( *alp );
    formula = std::string(Form("TMath::Max(1e-50,pow(%s,@%d))",polFormula.c_str(),argset.getSize()-1));

  }

  else if( funcFamily=="invpow" ) {

    RooRealVar* alp = new RooRealVar(Form("invpow%d_alpha" , iOrder), "invpow_alpha" , -4., -20, 0. );
    argset.add( *alp );
    formula = std::string(Form("TMath::Max(1e-50,pow(1+%s,@%d))",polFormula.c_str(),argset.getSize()-1));

  }

  else if( funcFamily=="invpowlin" ) {

    RooRealVar* c = new RooRealVar(Form("invpowlin%d_c" , iOrder), "invpowlin_c" , 0.001);
    c->setConstant(false);
    argset.add( *c );
    RooRealVar* alp = new RooRealVar(Form("invpowlin%d_alpha" , iOrder), "invpowlin_alpha" , -4.);
    alp->setConstant(false);
    argset.add( *alp );
    formula = std::string(Form("TMath::Max(1e-50,pow(1+@%d*@0, @%d+%s))", argset.getSize()-2, argset.getSize()-1, polFormula.c_str()));

  }

  else if( funcFamily=="moddijet" ) {

    RooRealVar* a = new RooRealVar(Form("moddijet%d_a" , iOrder), "moddijet_a" , 5., -100., 10.);
    a->setConstant(false);
    argset.add( *a );
    RooRealVar* b = new RooRealVar(Form("moddijet%d_b" , iOrder), "moddijet_b" , -1., -100., 10.);
    b->setConstant(false);
    argset.add( *b );
    RooRealVar* alp = new RooRealVar(Form("moddijet%d_alpha" , iOrder), "moddijet_alpha" , -50, -1000., 0.);
    alp->setConstant(false);
    argset.add( *alp );
    formula = std::string(Form("TMath::Max(1e-50,pow(@0,@%d+@%d*log(@0))*pow(1.-(%s),@%d))", argset.getSize()-3, argset.getSize()-2, polFormula.c_str(), argset.getSize()-1));

  }


  std::cout << formula << std::endl;
  RooGenericPdf* model = new RooGenericPdf( Form("%s_%d", funcFamily.c_str(), iOrder), funcFamily.c_str(), formula.c_str(), argset );

  RooFitResult *fitRes = model->fitTo(*data,Strategy(2),Warnings(false),Minimizer("Minuit2"),Save(true));
  fitRes->Print();
  thisNll = fitRes->minNll();

  RooPlot* frame = x->frame();
  data->plotOn(frame);
  model->plotOn(frame);
  frame->GetYaxis()->SetRangeUser( 0.01, 10000 );

  TCanvas* c1 = new TCanvas("c1", "c1", 600, 600 );
  c1->cd();
  c1->SetLogy();
  frame->Draw();

  c1->SaveAs(Form("%s/ftest_%s_o%d.eps", outdir.c_str(), funcFamily.c_str(),iOrder));
  c1->SaveAs(Form("%s/ftest_%s_o%d.pdf", outdir.c_str(), funcFamily.c_str(),iOrder));

  delete c1;
  
  return thisNll;

}




float fitWithModelSum( const std::string& outdir, RooRealVar* x, RooDataSet* data, const std::string& funcFamily, unsigned iOrder ) {

  float init = -4.;
  if( funcFamily=="exp" )
    init = -0.01;

  RooArgSet argset;
  argset.add( *x );
  std::string formula;
  getSumFunctionsFormula( funcFamily, iOrder, formula, argset, init );
  std::cout << formula << std::endl;

  RooGenericPdf* model = new RooGenericPdf( Form("TMath::Max(1e-50,%s_%d)", funcFamily.c_str(), iOrder), funcFamily.c_str(), formula.c_str(), argset );

  RooFitResult *fitRes = model->fitTo(*data,Strategy(2),Warnings(false),Minimizer("Minuit2"),Save(true));
  fitRes->Print();
  float thisNll = fitRes->minNll();

  RooPlot* frame = x->frame();
  data->plotOn(frame);
  model->plotOn(frame);
  frame->GetYaxis()->SetRangeUser( 0.01, 10000 );

  TCanvas* c1 = new TCanvas("c1", "c1", 600, 600 );
  c1->cd();
  c1->SetLogy();
  frame->Draw();

  c1->SaveAs(Form("%s/ftest_sum%s_o%d.eps", outdir.c_str(), funcFamily.c_str(),iOrder));
  c1->SaveAs(Form("%s/ftest_sum%s_o%d.pdf", outdir.c_str(), funcFamily.c_str(),iOrder));

  delete c1;
  
  return thisNll;

}


void getPol( const std::string& name, unsigned iOrder, std::string& polFormula, RooArgSet& polargset, float init, float xmin, float xmax ) {

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

      RooRealVar* par = new RooRealVar( Form("%s_o%d_p%d", name.c_str(), iOrder, index), Form("p%d", index), pow(init,index) );
      par->setConstant(false);
      if( xmin>-9999. && xmax>-9999. )
        par->setRange(xmin,xmax);
      polargset.add( *par );

    }

  } // for iorder

  std::cout << polFormula << std::endl;
  polargset.Print();

}


void getSumFunctionsFormula( const std::string& funcFamily, unsigned iOrder, std::string& formula, RooArgSet& argset, float init ) {

  formula = "";

  unsigned additionalPars = iOrder*2 - 1;

  for( unsigned i=1; i<=additionalPars; ++i ) {

    if( i!=1 )  {
      formula += std::string(Form("+@%d*", i));
      RooRealVar* par = new RooRealVar( Form("%s_o%d_p%d", funcFamily.c_str(), iOrder, i), Form("p%d", i), 1., 0., 100.);
      par->setConstant(false);
      argset.add( *par );
      i++;
    }

    if( funcFamily=="pow" )
      formula += std::string(Form("%s(@0,@%d)",funcFamily.c_str(),i));
    else
      formula += std::string(Form("%s(@0*@%d)",funcFamily.c_str(),i));
    RooRealVar* par = new RooRealVar( Form("%s_o%d_p%d", funcFamily.c_str(), iOrder, i), Form("p%d", i), init );
    par->setConstant(false);
    argset.add( *par );

  } // for i

}
