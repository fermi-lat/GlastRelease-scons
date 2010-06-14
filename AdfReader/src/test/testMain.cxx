#include <iostream>
#include "AdfEvent/AdfEvent.h"
#include "AncillaryDataEvent/Digi.h"

#include "AncillaryDataUtil/AncillaryDataServer.h"

// do something useful another day
char *fileName="../data/Ped_DAQBARI_114_330000723.bin";

AncillaryDataServer *m_dataServer;
AncillaryData::AdfEvent *m_adfEvent;
//AncillaryData::Digi *m_digiEvent;


int initialize()
{
  std::cout << "---> initialize()" << std::endl << std::endl;
  m_dataServer = new AncillaryDataServer(fileName);
  m_dataServer->open();
  m_dataServer->printInformation();
  return 1;
}  


int execute()
{
  std::cout << "---> execute()" << std::endl << std::endl;
  m_adfEvent = m_dataServer->getNextEvent();
  if(m_adfEvent==0) return 0;
  m_adfEvent->print();
  // Creates Digis...
  AncillaryData::Digi digiEvent(m_adfEvent);

  //  m_adfEvent->dump();
  return 1;
}

int finalize()
{
  std::cout << "---> finalize()" << std::endl << std::endl;
  m_dataServer->close();
  delete m_dataServer;
}



int main(){
  int sc=initialize();
  int ev=0;

  while(ev<10000 && sc) 
    {
      std::cout<<"Event # "<<ev++<<std::endl;
      sc=execute();
    }
  finalize();

  /*
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
  */
  return 0;
}


