#include <iostream>
#include "AdfEvent/AdfEvent.h"

// do something useful another day
int main(){

  AncillaryData::EventSummaryData d;
  std::cout<<"qui"<<std::endl;
  unsigned int word[4] = {(200<<20),0,0,0};
  d.setData(word);
  std::cout<<"qui"<<std::endl;
  AncillaryData::AdfEvent a;
  std::cout<<"qui"<<std::endl;
  a.setEventSummaryData(d);
  std::cout<<"qui"<<std::endl;
  std::cout<<a.getEventSummaryData().getVersion()<<std::endl;
  return 0;
}

