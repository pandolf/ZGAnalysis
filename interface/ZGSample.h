#ifndef ZGSample_h
#define ZGSample_h

#include <string>
#include <vector>




class ZGSample {


 public:

  ZGSample();
  ~ZGSample();


  static std::vector<ZGSample> loadSamples(const std::string& filename, const std::string& filter="", int idMin=-1, int idMax=-1);
  static std::vector<ZGSample> loadSamples(const std::string& filename, int idMin, int idMax=-1);

  // publica data members:
  std::string name;
  std::string sname;
  std::string dir;
  std::string file;
  int id;
  int nevents;
  float xsection;
  float lumi;
  float kfact;
  float filter;
  float scale1fb;
  float PU_avg_weight;

 private:


};

#endif
