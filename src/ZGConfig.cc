#include "../interface/ZGConfig.h"

#include <iostream>
#include <fstream>
#include <cmath>

#include "TH1F.h"



ZGConfig::ZGConfig( const std::string& name ) {

  name_ = name;

  std::string configFileName = "cfgs/" + name + ".txt";

  std::cout << std::endl;
  std::cout << "-> Reading config file: " << configFileName << std::endl;
  std::cout << std::endl;

  lumi_ = 0.;
  lumi_JetHT_ = 0.;
  lumi_HTMHT_ = 0.;
  lumi_SinglePhoton_ = 0.;
  lumi_DoubleEG_ = 0.;
  lumi_DoubleMu_ = 0.;

  mcSamples_ = "";
  sigSamples_ = "";
  dataSamples_ = "";
  additionalStuff_ = "";
  analysisType_ = "zg";
  selection_ = "v0";
  smearing_ = true;


  std::ifstream IN(configFileName.c_str());
  char buffer[200];
  char StringValue[1000];


  while( IN.getline(buffer, 200, '\n') ) {

    if (buffer[0] == '#') {
      continue; // Skip lines commented with '#'                                                                                                                                                                                 
    }

    std::cout << buffer << std::endl;

    char this_name_c[200];
    sscanf(buffer, "%s %s", this_name_c, StringValue);
    std::string this_name(this_name_c);

    if( this_name=="lumi" )
      lumi_ = atof(StringValue);
    else if( this_name=="lumi_JetHT" )
      lumi_JetHT_ = atof(StringValue);
    else if( this_name=="lumi_HTMHT" )
      lumi_HTMHT_ = atof(StringValue);
    else if( this_name=="lumi_SinglePhoton" )
      lumi_SinglePhoton_ = atof(StringValue);
    else if( this_name=="lumi_DoubleEG" )
      lumi_DoubleEG_ = atof(StringValue);
    else if( this_name=="lumi_DoubleMu" )
      lumi_DoubleMu_ = atof(StringValue);
    else if( this_name=="mcSamples" )
      mcSamples_ = std::string(StringValue);
    else if( this_name=="sigSamples" )
      sigSamples_ = std::string(StringValue);
    else if( this_name=="dataSamples" )
      dataSamples_ = std::string(StringValue);
    else if( this_name=="additionalStuff" )
      additionalStuff_ = std::string(StringValue);
    else if( this_name=="analysisType" )
      analysisType_ = std::string(StringValue);
    else if( this_name=="selection" )
      selection_ = std::string(StringValue);
    else if( this_name=="smearing" ) {
      std::string smearing_str = std::string(StringValue);
      if( smearing_str=="true"  ) smearing_ = true;
      if( smearing_str=="false" ) smearing_ = false;
    }

  } // while getline

  std::cout << std::endl;


  if( this->useMC() && lumi_ <=0. ) {
    std::cout << "[ZGConfig] ERROR!! If you process MC files you need to set a valid lumi value in your cfg!" << std::endl;
    std::cout << "Exiting..." << std::endl;
    exit(761);
  }


     
}


float ZGConfig::lumi_JetHT() const { 

  return this->defaultLumi(lumi_JetHT_);

}


float ZGConfig::lumi_HTMHT() const { 

  return this->defaultLumi(lumi_HTMHT_);

}


float ZGConfig::lumi_SinglePhoton() const { 

  return this->defaultLumi(lumi_SinglePhoton_);

}


float ZGConfig::lumi_DoubleEG() const { 

  return this->defaultLumi(lumi_DoubleEG_);

}


float ZGConfig::lumi_DoubleMu() const { 

  return this->defaultLumi(lumi_DoubleMu_);

}


float ZGConfig::defaultLumi( float lumi ) const {

  float returnLumi = (lumi>0.) ? lumi : lumi_; // if not over-written, return (common) lumi_
  return returnLumi;

}



bool ZGConfig::useMC() const {

  return mcSamples_!="";

}



bool ZGConfig::dummyAnalysis() const {

  return dataSamples_=="datatest";

}


std::string ZGConfig::getEventYieldDir() const {

  std::string outputdir = "EventYields_" + name_;
  if( this->dummyAnalysis() ) outputdir += "_dummy";

  //double intpart;
  //double fracpart = modf(lumi_, &intpart);
  //std::string suffix;
  //if( fracpart>0. )
  //  suffix = std::string( Form("%.0fp%.0ffb", intpart, 10.*fracpart ) );
  //else
  //  suffix = std::string( Form("%.0ffb", intpart ) );
  //outputdir += suffix;

  return outputdir;

}




void ZGConfig::saveAs( const std::string& filename ) const {


  std::ofstream ofs(filename.c_str());

  ofs << "#name " << name_ << std::endl;

  ofs << "lumi "  << lumi_  << std::endl;
  ofs << "lumi_JetHT "        <<   lumi_JetHT_          << std::endl;
  ofs << "lumi_HTMHT "        <<   lumi_HTMHT_          << std::endl;
  ofs << "lumi_SinglePhoton " <<   lumi_SinglePhoton_   << std::endl;
  ofs << "lumi_DoubleEG "     <<   lumi_DoubleEG_       << std::endl;
  ofs << "lumi_DoubleMu "     <<   lumi_DoubleMu_       << std::endl;

  if( mcSamples_!="" )       ofs << "mcSamples " << mcSamples_ << std::endl;
  if( sigSamples_!="" )      ofs << "sigSamples " << sigSamples_ << std::endl;
  if( dataSamples_!="" )     ofs << "dataSamples " << dataSamples_ << std::endl;
  if( additionalStuff_!="" ) ofs << "additionalStuff " << additionalStuff_ << std::endl;
  if( selection_!="" ) ofs << "selection " << selection_ << std::endl;
  ofs << "smearing " << smearing_ << std::endl;


  std::cout << "[ZGConfig] Saved config file as '" << filename << "'." << std::endl;

}
