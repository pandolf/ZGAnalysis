#include <iostream>
#include <fstream>
#include <set>

using namespace std;

class EventKey {
public:
  EventKey(Int_t input_run=0, Int_t input_lumi=0, UInt_t input_evt=0) : 
    run_(input_run), lumi_(input_lumi), evt_(input_evt){;}

  Int_t run() const {return run_;}
  Int_t lumi() const {return lumi_;}
  UInt_t evt() const {return evt_;}

  bool operator<(EventKey const& right) const{
    if (run_ == right.run()) {
	return evt_ < right.evt();
    }
    return run_ < right.run();
  }

//  bool operator<(EventKey const& right) const{
//    if (run_ == right.run()) {
//      if (lumi_ == right.lumi()) {
//	return evt_ < right.evt();
//      }
//      return lumi_ < right.lumi();
//    }
//    return run_ < right.run();
//  }

private:
  Int_t run_;
  Int_t lumi_;
  UInt_t evt_;

};



void checkOverlap(string eventList_eth="events_eth.txt",
		      string eventList_cnc="events_cnc.txt") {

  std::set<EventKey> eth_events;
  ifstream ifs_eth(eventList_eth);
  std::cout << "Opening: " << eventList_eth << std::endl;
  while( ifs_eth.good() ) {
    int run, lumi;
    UInt_t evt;
    ifs_eth >> run >> lumi >> evt;
    EventKey newEvent(run,lumi,evt);
    eth_events.insert(newEvent);
  }
  std::cout << "Total ETH events: " << eth_events.size() << std::endl;


  std::set<EventKey> cnc_events;
  std::set<EventKey> overlap_events;
  ifstream ifs_cnc(eventList_cnc);
  std::cout << "Opening: " << eventList_cnc << std::endl;
  while( ifs_cnc.good() ) {
    int run, lumi;
    UInt_t evt;
    ifs_cnc >> run >> lumi >> evt;
    EventKey newEvent(run,lumi,evt);
    cnc_events.insert(newEvent);
    if( eth_events.find(newEvent)!=eth_events.end() )
      overlap_events.insert( newEvent );
  }

  std::cout << "Total CNC events: " << cnc_events.size() << std::endl;
  std::cout << "Overlap between the two analyses: " << overlap_events.size() << "/" << eth_events.size() << " (" << (float)overlap_events.size()/((float)eth_events.size())*100. << "%%)" << std::endl;

}
