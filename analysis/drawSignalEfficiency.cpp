#include <iostream>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TLegend.h"

#include "../interface/ZGConfig.h"
#include "../interface/ZGSample.h"
#include "../interface/ZGDrawTools.h"



void addEfficiencyPoint( TGraphErrors* gr_eff, const ZGSample& sample, const ZGConfig& cfg, const std::string& sel="" );
TGraphErrors* getRatio( TGraphErrors* gr, TF1* f1 );
void drawCompare( const ZGConfig& cfg, const std::string& width, TGraphErrors* gr_eff_ee, TGraphErrors* gr_eff_mm, TF1* f1_ee, TF1* f1_mm );
void drawCompare( const ZGConfig& cfg, const std::string& width, TGraphErrors* gr_eff, TF1* f1, const std::string& name, const std::string& channel );
TF1* drawRelativeFractions( const ZGConfig& cfg, const std::string& width, TGraphErrors* gr_eff_ee, TGraphErrors* gr_eff_mm );
TF1* fitConstantPar( TF1* f1, TGraphErrors* gr_eff, float shift=-999. );
TF1* fitEfficiency( TGraphErrors* gr );
TString widthText( const std::string& width );



float xMax = 2100.;

int main( int argc, char* argv[] ) {

  if( argc<2 ) {
    std::cout << "USAGE: ./drawSignalEfficiency [configFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  ZGDrawTools::setStyle();

  std::string configFileName(argv[1]);
  ZGConfig cfg(configFileName);

  std::string samplesFileName = "../samples/samples_" + cfg.mcSamples() + ".dat";
  std::cout << std::endl << std::endl;
  std::cout << "-> Loading signal samples from file: " << samplesFileName << std::endl;

  std::vector<ZGSample> fSamples = ZGSample::loadSamples(samplesFileName, 1000, 10000); // only signal (id>=1000)


  if( fSamples.size()==0 ) {
    std::cout << "No signal samples found." << std::endl;
    exit(111);
  } 


  std::vector<std::string> widths;
  widths.push_back("0p014");
  widths.push_back("5p6");

  for( unsigned iWidth=0;iWidth<widths.size(); ++iWidth ) {

    //if( widths[iWidth]=="5p6" ) xMax = 1000.;
    //else                        xMax = 2100.;

    TFile* file = TFile::Open( Form("%s/signalEfficiency_w%s.root", cfg.getEventYieldDir().c_str(), widths[iWidth].c_str()), "recreate" );
    file->cd();

    TGraphErrors* gr_eff    = new TGraphErrors(0);
    TGraphErrors* gr_eff_ee = new TGraphErrors(0);
    TGraphErrors* gr_eff_mm = new TGraphErrors(0);
    gr_eff   ->SetName("eff_all");
    gr_eff_ee->SetName("eff_ee");
    gr_eff_mm->SetName("eff_mm");

    gr_eff_ee->SetTitle("e^{+}e^{-}#gamma");
    gr_eff_mm->SetTitle("#mu^{+}#mu^{-}#gamma");

    for( unsigned i=0; i<fSamples.size(); ++i ) {
      TString name_tstr(fSamples[i].name);
      if( !name_tstr.Contains(widths[iWidth]) ) continue;
      std::cout << "-> Starting sample: " << fSamples[i].name << std::endl;
      addEfficiencyPoint( gr_eff   , fSamples[i], cfg );
      addEfficiencyPoint( gr_eff_ee, fSamples[i], cfg, "ee" );
      addEfficiencyPoint( gr_eff_mm, fSamples[i], cfg, "mm" );
    }
    


    //// compare to eff x acc:
    //TFile* file_aXe = TFile::Open("genAcceptanceTimesEfficiency.root");
    //TF1* f1_aXe_0p014    = (TF1*)file_aXe->Get("f1_gr_0p014_times_line_all");
    //TF1* f1_aXe_0p014_ee = (TF1*)file_aXe->Get("f1_gr_0p014_times_line_ee");
    //TF1* f1_aXe_0p014_mm = (TF1*)file_aXe->Get("f1_gr_0p014_times_line_mm");

    //drawCompare( cfg, gr_eff   , f1_aXe_0p014   , "all", "Electrons + Muons" );
    //drawCompare( cfg, gr_eff_ee, f1_aXe_0p014_ee, "ee", "Electron Channel" );
    //drawCompare( cfg, gr_eff_mm, f1_aXe_0p014_mm, "mm", "Muon Channel" );

    ////TF1* f1_aXe_0p014_scale    = fitConstantPar( f1_aXe_0p014   , gr_eff   , 0.0285738 );
    ////TF1* f1_aXe_0p014_ee_scale = fitConstantPar( f1_aXe_0p014_ee, gr_eff_ee, 0.04 );
    ////TF1* f1_aXe_0p014_mm_scale = fitConstantPar( f1_aXe_0p014_mm, gr_eff_mm, 0.0 );
    //TF1* f1_aXe_0p014_scale    = fitConstantPar( f1_aXe_0p014   , gr_eff    );
    //TF1* f1_aXe_0p014_ee_scale = fitConstantPar( f1_aXe_0p014_ee, gr_eff_ee );
    //TF1* f1_aXe_0p014_mm_scale = fitConstantPar( f1_aXe_0p014_mm, gr_eff_mm );

    //drawCompare( cfg, gr_eff   , f1_aXe_0p014_scale   , "all_scale", "Electrons + Muons" );
    //drawCompare( cfg, gr_eff_ee, f1_aXe_0p014_ee_scale, "ee_scale", "Electron Channel" );
    //drawCompare( cfg, gr_eff_mm, f1_aXe_0p014_mm_scale, "mm_scale", "Muon Channel" );

    TF1* f1_aXe_0p014_cust    = fitEfficiency( gr_eff    );
    TF1* f1_aXe_0p014_ee_cust = fitEfficiency( gr_eff_ee );
    TF1* f1_aXe_0p014_mm_cust = fitEfficiency( gr_eff_mm );

    drawCompare( cfg, widths[iWidth], gr_eff   , f1_aXe_0p014_cust   , "all", std::string(widthText(widths[iWidth]) + " (ee#gamma+#mu#mu#gamma)") );
    drawCompare( cfg, widths[iWidth], gr_eff_ee, f1_aXe_0p014_ee_cust, "ee" , std::string(widthText(widths[iWidth]) + " (ee#gamma)") );
    drawCompare( cfg, widths[iWidth], gr_eff_mm, f1_aXe_0p014_mm_cust, "mm" , std::string(widthText(widths[iWidth]) + " (#mu#mu#gamma)") );

    // compare electrons and muons
    drawCompare( cfg, widths[iWidth], gr_eff_ee, gr_eff_mm, f1_aXe_0p014_ee_cust, f1_aXe_0p014_mm_cust );
    TF1* f1_frac_ee = drawRelativeFractions( cfg, widths[iWidth], gr_eff_ee, gr_eff_mm );


    file->cd();
    gr_eff->Write();
    gr_eff_ee->Write();
    gr_eff_mm->Write();
    //f1_aXe_0p014_scale   ->SetName("f1_0p014_all");
    //f1_aXe_0p014_ee_scale->SetName("f1_0p014_ee");
    //f1_aXe_0p014_mm_scale->SetName("f1_0p014_mm");
    //f1_aXe_0p014_scale   ->Write();
    //f1_aXe_0p014_ee_scale->Write();
    //f1_aXe_0p014_mm_scale->Write();
    f1_aXe_0p014_cust   ->SetName("f1_eff_all");
    f1_aXe_0p014_ee_cust->SetName("f1_eff_ee");
    f1_aXe_0p014_mm_cust->SetName("f1_eff_mm");
    f1_aXe_0p014_cust   ->Write();
    f1_aXe_0p014_ee_cust->Write();
    f1_aXe_0p014_mm_cust->Write();
    f1_frac_ee->SetName("f1_frac_ee");
    f1_frac_ee->Write();
    file->Close();

    std::cout << "-> Wrote signal efficiencies for width = " << widths[iWidth] << " in: " << file->GetName() << std::endl;

  }


  return 0;

}



void drawCompare( const ZGConfig& cfg, const std::string& width, TGraphErrors* gr_eff_ee, TGraphErrors* gr_eff_mm, TF1* f1_ee, TF1* f1_mm ) {

  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();
  
  TH2D* h2_axes = new TH2D("axes", "", 10, 300., xMax, 10, 0., 1.0001 );
  h2_axes->SetXTitle( "Generated Z#gamma Mass [GeV]" );
  h2_axes->GetXaxis()->SetNdivisions(1006, false);
  h2_axes->SetYTitle( "Efficiency #times Acceptance" );
  h2_axes->Draw();

  gr_eff_ee->SetMarkerStyle(20);
  gr_eff_ee->SetMarkerSize(1.3);
  gr_eff_ee->SetMarkerColor(46);
  gr_eff_ee->SetLineColor(46);
  gr_eff_ee->SetLineWidth(2);

  gr_eff_mm->SetMarkerStyle(20);
  gr_eff_mm->SetMarkerSize(1.3);
  gr_eff_mm->SetMarkerColor(38);
  gr_eff_mm->SetLineColor(38);
  gr_eff_mm->SetLineWidth(2);

  if( f1_ee!=0 ) {
    f1_ee->SetLineColor(46);
    f1_ee->Draw("same");
  }
  if( f1_mm!=0 ) {
    f1_mm->SetLineColor(38);
    f1_mm->Draw("same");
  }

  gr_eff_ee->Draw("p same");
  gr_eff_mm->Draw("p same");

  TLegend* legend = new TLegend( 0.2, 0.7, 0.55, 0.9 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  legend->SetHeader( widthText(width) );
  legend->AddEntry( gr_eff_mm, gr_eff_mm->GetTitle(), "PL" );
  legend->AddEntry( gr_eff_ee, gr_eff_ee->GetTitle(), "PL" );
  legend->Draw("same");

  ZGDrawTools::addLabels( c1, -1, "CMS Simulation" );

  c1->SaveAs( Form("%s/signalEfficiency_w%s_%s_vs_%s.eps", cfg.getEventYieldDir().c_str(), width.c_str(), gr_eff_ee->GetName(), gr_eff_mm->GetName()) );
  c1->SaveAs( Form("%s/signalEfficiency_w%s_%s_vs_%s.pdf", cfg.getEventYieldDir().c_str(), width.c_str(), gr_eff_ee->GetName(), gr_eff_mm->GetName()) );

  c1->SetLogx();

  c1->SaveAs( Form("%s/signalEfficiency_w%s_%s_vs_%s_logx.eps", cfg.getEventYieldDir().c_str(), width.c_str(), gr_eff_ee->GetName(), gr_eff_mm->GetName()) );
  c1->SaveAs( Form("%s/signalEfficiency_w%s_%s_vs_%s_logx.pdf", cfg.getEventYieldDir().c_str(), width.c_str(), gr_eff_ee->GetName(), gr_eff_mm->GetName()) );

  delete legend;
  delete c1;
  delete h2_axes;

}




TF1* drawRelativeFractions( const ZGConfig& cfg, const std::string& width, TGraphErrors* gr_eff_ee, TGraphErrors* gr_eff_mm ) {


  TGraphErrors* gr_ee_rel = new TGraphErrors(0);
  TGraphErrors* gr_mm_rel = new TGraphErrors(0);

  for( int iPoint=0; iPoint<gr_eff_ee->GetN(); ++iPoint ) {

    Double_t x_ee, eff_ee;
    gr_eff_ee->GetPoint( iPoint, x_ee, eff_ee );

    Double_t x_mm, eff_mm;
    gr_eff_mm->GetPoint( iPoint, x_mm, eff_mm );

    Double_t eff_ee_err = gr_eff_ee->GetErrorY( iPoint );
    Double_t eff_mm_err = gr_eff_mm->GetErrorY( iPoint );

    float tot = eff_ee + eff_mm;
    gr_ee_rel->SetPoint( iPoint, x_ee, eff_ee/tot );
    gr_mm_rel->SetPoint( iPoint, x_mm, eff_mm/tot );
    gr_ee_rel->SetPointError( iPoint, 0., eff_ee_err );
    gr_mm_rel->SetPointError( iPoint, 0., eff_mm_err );

  }

  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();
  
  TH2D* h2_axes = new TH2D("axes", "", 10, 300., xMax, 10, 0., 1.0001 );
  h2_axes->SetXTitle( "Generated Z#gamma Mass [GeV]" );
  h2_axes->SetYTitle( "Relative Amount" );
  h2_axes->Draw();


  TF1* line_ee = new TF1("line_ee", "[0]+[1]*x", 300., xMax);
  line_ee->SetLineColor(46);
  gr_ee_rel->Fit( line_ee, "QR" );

  TF1* line_mm = new TF1("line_mm", "[0]+[1]*x", 300., xMax);
  line_mm->SetLineColor(38);
  gr_mm_rel->Fit( line_mm, "QR" );

  gr_ee_rel->SetMarkerStyle(20);
  gr_ee_rel->SetMarkerColor(46);
  gr_ee_rel->SetMarkerSize(1.3);
  gr_ee_rel->SetLineColor(46);
  gr_ee_rel->SetLineWidth(2);

  gr_mm_rel->SetMarkerStyle(20);
  gr_mm_rel->SetMarkerColor(38);
  gr_mm_rel->SetMarkerSize(1.3);
  gr_mm_rel->SetLineColor(38);
  gr_mm_rel->SetLineWidth(2);

  gr_mm_rel->Draw("p same");
  gr_ee_rel->Draw("p same");

  TLegend* legend = new TLegend( 0.2, 0.7, 0.55, 0.9 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  //legend->AddEntry( gr_mm_rel, Form("#mu^{+}#mu^{-}#gamma (%.1f %%)", line_mm->GetParameter(0)*100.), "PL" );
  //legend->AddEntry( gr_ee_rel, Form("e^{+}e^{-}#gamma (%.1f %%)"    , line_ee->GetParameter(0)*100.), "PL" );
  legend->AddEntry( gr_mm_rel, "#mu^{+}#mu^{-}#gamma", "PL" );
  legend->AddEntry( gr_ee_rel, "e^{+}e^{-}#gamma"    , "PL" );
  legend->Draw("same");

  ZGDrawTools::addLabels( c1, -1, "CMS Simulation" );

  c1->SaveAs( Form("%s/signalEfficiency_rel_w%s_ee_vs_mm.eps", cfg.getEventYieldDir().c_str(), width.c_str()) );
  c1->SaveAs( Form("%s/signalEfficiency_rel_w%s_ee_vs_mm.pdf", cfg.getEventYieldDir().c_str(), width.c_str()) );


  delete legend;
  delete c1;
  delete h2_axes;
  delete line_mm;
  delete gr_ee_rel;
  delete gr_mm_rel;

  return line_ee;

}





void drawCompare( const ZGConfig& cfg, const std::string& width, TGraphErrors* gr_eff, TF1* f1, const std::string& name, const std::string& channel ) {

  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();
  
  TH2D* h2_axes = new TH2D("axes", "", 10, 300., xMax, 10, 0., 1.0001 );
  h2_axes->SetXTitle( "Generated Z#gamma Mass [GeV]" );
  h2_axes->SetYTitle( "Efficiency #times Acceptance" );
  h2_axes->Draw();

  f1->SetLineColor(46);
  f1->Draw("same");

  gr_eff->SetMarkerStyle(20);
  gr_eff->SetMarkerColor(kBlack);
  gr_eff->SetMarkerSize(1.3);
  gr_eff->Draw("p same");

  TLegend* legend = new TLegend( 0.2, 0.7, 0.55, 0.9, channel.c_str() );
  legend->SetFillColor(0);
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);
  legend->AddEntry( gr_eff, "From Full Sim Samples", "P" );
  legend->AddEntry( f1, "Fit", "L" );
  legend->Draw("same");

  ZGDrawTools::addLabels( c1, -1, "CMS Simulation" );

  c1->SaveAs( Form("%s/signalEfficiency_w%s_%s.eps", cfg.getEventYieldDir().c_str(), width.c_str(), name.c_str()) );
  c1->SaveAs( Form("%s/signalEfficiency_w%s_%s.pdf", cfg.getEventYieldDir().c_str(), width.c_str(), name.c_str()) );


  //TGraphErrors* gr_ratio = getRatio( gr_eff, f1);
  //TF1* line_ratio = new TF1("line_ratio", "[0]+[1]*x", 400., 1000.);
  //gr_ratio->Fit( line_ratio, "QR" );

  //float ratio = line_ratio->GetParameter(0);

  //delete line_ratio;
  //delete gr_ratio;
  delete legend;
  delete c1;
  delete h2_axes;

  //return ratio;

}



void addEfficiencyPoint( TGraphErrors* gr_eff, const ZGSample& sample, const ZGConfig& cfg, const std::string& sel ) {

  std::string name = sample.name;
  TString name_tstr(name);

  float nTotalGenEvents = (float)sample.nevents;
  nTotalGenEvents = nTotalGenEvents*2./3.; // fuck taus
  if( sel=="ee" || sel=="mm" ) nTotalGenEvents/=2.;

  TFile* file_pass = TFile::Open( Form("%s/trees.root", cfg.getEventYieldDir().c_str()) );
  TTree* tree_pass = (TTree*)file_pass->Get( name.c_str() );
  int nPassEvents;
  if( sel!="" ) {
    std::string tree_sel="";
    if( sel=="ee" ) tree_sel = "leptType==11";
    if( sel=="mm" ) tree_sel = "leptType==13";
    nPassEvents =  tree_pass->GetEntries(Form("%s", tree_sel.c_str()));
    //nPassEvents =  tree_pass->GetEntries(Form("passHLT && %s", tree_sel.c_str()));
  } else {
    nPassEvents =  tree_pass->GetEntries("");
    //nPassEvents =  tree_pass->GetEntries("passHLT");
  }

  float thisEff = (float)nPassEvents/((float)nTotalGenEvents);
  float thisEffErr = sqrt( thisEff*(1.-thisEff)/((float)nTotalGenEvents) );

  int mass = 0;
  if( name_tstr.BeginsWith("GluGlu") ) {
    std::string delimiter = "M";
    std::string mass_str = name.substr(name.find(delimiter)+delimiter.length(), 4);
    mass = atoi(mass_str.c_str());
  } else {
    std::string delimiter = "M_";
    std::string mass_str = name.substr(name.find(delimiter)+delimiter.length(), 4);
    mass = atoi(mass_str.c_str());
  }

  if( mass > xMax || mass==0 ) return;

  int iPoint = gr_eff->GetN();
  gr_eff->SetPoint( iPoint, mass, thisEff );
  gr_eff->SetPointError( iPoint, 0., thisEffErr );

}


TGraphErrors* getRatio( TGraphErrors* gr, TF1* f1 ) {


  TGraphErrors* gr_ratio = new TGraphErrors(0);

  for( int iPoint=0; iPoint<gr->GetN(); ++iPoint ) {

    Double_t x, y;
    gr->GetPoint( iPoint, x, y );
    Double_t y_err = gr->GetErrorY( iPoint );

    gr_ratio->SetPoint( iPoint, x, y/f1->Eval(x) );
    gr_ratio->SetPointError( iPoint, 0., y_err/f1->Eval(x) );

  }

  return gr_ratio;

}



TF1* fitConstantPar( TF1* f1, TGraphErrors* gr_eff, float shift ) {

  Double_t xmin, xmax;
  f1->GetRange( xmin, xmax );
  TF1* f1_new = new TF1( Form("%s_scale", f1->GetName()), f1->GetExpFormula(), xmin, xmax );

  for( int ipar=0; ipar<f1->GetNpar(); ++ipar ) {
    if( ipar!=0 )
      f1_new->FixParameter( ipar, f1->GetParameter(ipar) );
    else
      f1_new->SetParameter( ipar, f1->GetParameter(ipar) );
  }

  if( shift>-999. ) {
    f1_new->SetParameter(0, f1_new->GetParameter(0) + shift);
  } else {
    gr_eff->Fit( f1_new, "QR0" );
  }

  float diff = f1_new->GetParameter(0)-f1->GetParameter(0);
  std::cout << "Added " << diff << " to: " << f1->GetName() << std::endl;

  return f1_new;

}


TF1* fitEfficiency( TGraphErrors* gr ) {

  TF1* f1_new = new TF1( Form("f1_%s_cust", gr->GetName()), "[0] + [1]*(x-2000.)*(x-2000.) + [2]*(x-2000.)*(x-2000.)*(x-2000.)*(x-2000.)", 300., xMax );
  //TF1* f1_new = new TF1( Form("f1_%s_cust", gr->GetName()), "[0] + [1]*(x-2000.)*(x-2000.)", 300., xMax );
  //TF1* f1_new = new TF1( Form("f1_%s_cust", gr->GetName()), "[0] + [1]*x + [2]*x*x", 300., xMax );

  gr->Fit( f1_new, "QR0" );

  return f1_new;

}


TString widthText( const std::string& width ) {

  TString width_tstr(width);
  width_tstr.ReplaceAll( "p", "." );

  width_tstr = "W = " + width_tstr + "%";
  
  return width_tstr;

}
