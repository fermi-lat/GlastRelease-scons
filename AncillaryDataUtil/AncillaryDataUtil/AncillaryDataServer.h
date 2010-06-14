#ifndef ANCILLARYDATASERVER_HH
#define ANCILLARYDATASERVER_HH

#include "AncillaryDataUtil/AncillaryDataHeader.h"
#include "AncillaryDataUtil/AncillaryDataTail.h"
#include "AdfEvent/AdfEvent.h"
#include "AdfEvent/EventSummaryData.h"
using namespace AncillaryData;

class AncillaryDataServer
{
 public:
 
  AncillaryDataServer(std::string dataFileName);
  ~AncillaryDataServer();
  void open();
  void close();
  void synchronizationLost();
  int  readRawWord();
  AncillaryWord getNextDataWord();

  void readFileHeader();
  bool  readEventHeader(AdfEvent *currentEvent);
  bool  readEventData(AdfEvent *currentEvent);
  AdfEvent *getNextEvent();
  void printInformation();
  
  
 private:
 
  std::string m_dataFileName; 
  FILE *m_dataFile;
  AncillaryData::AncillaryDataHeader m_header;
  AncillaryData::AdfEvent            m_event;
  AncillaryData::AncillaryDataTail   m_tail;
};
#endif
