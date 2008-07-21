#include <iostream>

#include "facilities/Util.h"
#include "AdfEvent/AdfEvent.h"
#include "AncillaryDataEvent/Digi.h"
#include "AncillaryDataUtil/AncillaryDataServer.h"

// do something useful another day
//using namespace AncillaryData
int main()
{  
  std::string dataFilePath="$(ADFREADERDATAPATH)/CR_01_v3.bin";
  facilities::Util::expandEnvVar(&dataFilePath);
  AncillaryDataServer *m_dataServer = new AncillaryDataServer(dataFilePath);
  m_dataServer->open();
  m_dataServer->printInformation();
  for (int i =0; i<100; i++)
    {
      AncillaryData::AdfEvent *currentEvent = m_dataServer->getNextEvent();
      currentEvent->print();
      currentEvent->dump();
      //////////////////////////////////////////////////
      AncillaryData::Digi *digiEvent = new AncillaryData::Digi(currentEvent);
      digiEvent->print();
    }
  m_dataServer->close();   
  return 0;
}
