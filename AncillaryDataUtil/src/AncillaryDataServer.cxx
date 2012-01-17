#include "AdfEvent/AncillaryWord.h"
#include "AncillaryDataUtil/AncillaryDataTail.h"
#include "AncillaryDataUtil/AncillaryDataServer.h"
#include <cstdlib>
#include <cstdio>

#define DEBUG 0

using namespace AncillaryData;

AncillaryDataServer::AncillaryDataServer(std::string dataFileName)
{
  m_dataFileName = dataFileName;
  m_counter=0;
}

AncillaryDataServer::~AncillaryDataServer()
{
  //  if (m_currentEvent) delete m_currentEvent;
}

void AncillaryDataServer::open()
{
  if(!(m_dataFile = fopen(m_dataFileName.c_str(), "r"))){
    std::cout << "Binary data file not found." << std::endl;
    std::cout << "Exiting..." << std::endl;
    exit(0);
  }
  readFileHeader();
}

void AncillaryDataServer::close()
{
  fclose(m_dataFile);
}

void AncillaryDataServer::synchronizationLost()
{
  std::cout << "Syncrhonization lost.\n Exiting..." << std::endl;
  exit(0);
}

int AncillaryDataServer::readRawWord()
{
  int word[1];
  fread(word, ANCILLARY_WORD_LENGTH, 1, m_dataFile);
  if (DEBUG) std::cout<<"--> "<<std::hex<<word[0]<<std::endl;
  return word[0];
}


AncillaryWord AncillaryDataServer::getNextDataWord()
{
  AncillaryWord dataWord;
  dataWord.setData(readRawWord());
  return dataWord;
}

void AncillaryDataServer::readFileHeader()
{ 
  unsigned int word[HEADER_LENGTH];
  fread(word, ANCILLARY_WORD_LENGTH, HEADER_LENGTH, m_dataFile);
  m_header.setData(word);
  fseek(m_dataFile, m_header.getHeaderLength(), SEEK_SET);
}	 

bool  AncillaryDataServer::readEventHeader(AdfEvent *currentEvent)
{ 
  unsigned word[4];
  EventSummaryData summaryData;
  fread(word, ANCILLARY_WORD_LENGTH, 4, m_dataFile);
  summaryData.setData(word);
  currentEvent->setEventSummaryData(summaryData);    
  if(!summaryData.checkVersion()) 
    {
      std::cout<<" ******************** WARNING WRONG EVENT HEADER VERSION !! WE KEEP GOING..."<<std::endl;
      std::cout<<word[0]<<std::endl;
      std::cout<<word[1]<<std::endl;
      std::cout<<word[2]<<std::endl;
      std::cout<<word[3]<<std::dec<<std::endl;
      return false;
    }
  return true;
}

bool  AncillaryDataServer::readEventData(AdfEvent *currentEvent)
{ 
  //  AncillaryWord OneWord;
  unsigned remainingBytes = currentEvent->getEventSummaryData().getEventLength(); 
  remainingBytes-=4*ANCILLARY_WORD_LENGTH;
  if(remainingBytes<=0) return false;
  while(remainingBytes>0)
    {
      unsigned int rawword = readRawWord();
      remainingBytes-=ANCILLARY_WORD_LENGTH;
      AncillaryWord OneWord;
      OneWord.setData(rawword);
      //      AncillaryWord OneWord = getNextDataWord();
      unsigned int Header = OneWord.getHeader();
      QdcHeaderWord qdcHW;
      FadcHeaderWord fadcHW;
      ScalerHeaderWord scalerHW;
      
      if(DEBUG) std::cout<<"remainingBytes "<<remainingBytes<<std::endl;

      switch (Header)
	{
	case ANCILLARY_QDC_HID:
	  qdcHW.setData(rawword);
	  if(DEBUG) std::cout<<"Check QDC Header: "<<qdcHW.checkHeader()<<" "<<qdcHW.getQdcFifo()<<std::endl;

	  for(unsigned int i=0; i < qdcHW.getQdcFifo(); i++)
	    {
	      QdcDataWord qdcDW;
	      qdcDW.setData(readRawWord());
	      qdcDW.checkHeader();
	      currentEvent->appendQdcDataWord(qdcDW);
	    }
	  remainingBytes-=qdcHW.getQdcFifo()*ANCILLARY_WORD_LENGTH;
	  break;
	case ANCILLARY_FADC_HID:
	  fadcHW.setData(rawword);
	     if(DEBUG) std::cout<<"Check FADC Header: "<<fadcHW.checkHeader()<<" "<<fadcHW.getFadcFifo()<<std::endl;
	  for(unsigned int i=0; i < fadcHW.getFadcFifo(); i++)
	    {
	      FadcDataWord fadcDW;
	      fadcDW.setData(readRawWord());
	      fadcDW.checkHeader();
	      currentEvent->appendFadcDataWord(fadcDW);
	    }
	  remainingBytes-=fadcHW.getFadcFifo()*ANCILLARY_WORD_LENGTH;
	  break;
	case ANCILLARY_SCALER_ID:
	  scalerHW.setData(rawword);
	  if(DEBUG) std::cout<<"Check SCALER Header: "<<scalerHW.checkHeader()<<" "<<scalerHW.getScalerCounters()<<std::endl;
	  for(unsigned int i=0; i < scalerHW.getScalerCounters(); i++)
	    {
	      ScalerDataWord scalerDW;
	      scalerDW.setData(readRawWord());
	      currentEvent->appendScalerDataWord(scalerDW);
	    }
	  remainingBytes-= scalerHW.getScalerCounters()*ANCILLARY_WORD_LENGTH;
	  break;
	default:
	  break;
	}
    }
  return true;
}
  
AdfEvent *AncillaryDataServer::getNextEvent()
{

  if(m_counter%10==0 && DEBUG) std::cout<<" Get event number "<<++m_counter<<std::endl;
  // Get the event header:
  AdfEvent *currentEvent = new AdfEvent();
  bool suxess;
  suxess = readEventHeader(currentEvent);
  if(!suxess) 
    {
      std::cout<<"failed to read the Event header: "<<std::endl;
      suxess=readFileTail();
      if(!suxess) 
	{
	  std::cout<<"failed to read the TAIL "<<std::endl;
	}
      return 0;
    }
  suxess = readEventData(currentEvent);
  if (suxess) return currentEvent;
  std::cout<<"failed to read the Event Data: "<<std::endl;
  return 0;
}

void AncillaryDataServer::printInformation()
{
  std::cout << "#-----------------------------------------------"          << std::endl;
  std::cout << "Ancillary data server information:"                       << std::endl;
  std::cout << "Version:                " << m_header.getVersion()        << std::endl;
  std::cout << "Header Length (bytes):  " << m_header.getHeaderLength()   << std::endl;
  std::cout << "Run Number:             " << m_header.getRunNumber()      << std::endl;
  std::cout << "Comment Length (bytes): " << m_header.getCommentLength()  << std::endl;
  std::cout << "#-----------------------------------------------"          << std::endl;
}

bool AncillaryDataServer::readFileTail()
{
  unsigned int word[TAIL_LENGTH];
  int count = fread(word, ANCILLARY_WORD_LENGTH, TAIL_LENGTH, m_dataFile);
  if (count != TAIL_LENGTH) return(false); // HMK actually check to see all is well
  m_tail.setData(word);
  std::cout << "FileTail, Read: "<<m_tail.totalNumEvents()<<std::endl;
  return(true);  // Need to return a value for windows compile - assuming true HMK
}
