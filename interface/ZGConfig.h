#ifndef ZGConfig_h
#define ZGConfig_h

#include <string>

class ZGConfig {

 public:

  ZGConfig( const std::string& configFileName );

  std::string name() const { return name_; };

  float lumi()      const { return lumi_; };
  float lumi_JetHT()         const;
  float lumi_HTMHT()         const;
  float lumi_SinglePhoton()  const;
  float lumi_DoubleEG()      const;
  float lumi_DoubleMu()      const;

 
  std::string mcSamples()       const { return mcSamples_; };
  std::string sigSamples()      const { return sigSamples_; };
  std::string dataSamples()     const { return dataSamples_; };
  std::string additionalStuff() const { return additionalStuff_; };
  std::string analysisType()    const { return analysisType_; };
  std::string selection()       const { return selection_; };
  bool smearing()       const { return smearing_; };


  bool useMC() const;

  bool dummyAnalysis() const;

  std::string getEventYieldDir() const;

  void saveAs( const std::string& filename ) const;

 private:

  float defaultLumi( float lumi) const;

  std::string name_;

  float lumi_;
  float lumi_JetHT_;
  float lumi_HTMHT_;
  float lumi_SinglePhoton_;
  float lumi_DoubleEG_;
  float lumi_DoubleMu_;

  std::string mcSamples_;
  std::string sigSamples_;
  std::string dataSamples_;
  std::string additionalStuff_;
  std::string analysisType_;
  std::string selection_;
  bool smearing_;

};



#endif
