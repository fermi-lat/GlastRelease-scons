#include "AdfEvent/AncillaryWord.h"
#include "AncillaryDataUtil/AncillaryDataServer.h"

using namespace AncillaryData;

AncillaryDataServer::AncillaryDataServer(std::string dataFileName)
{
  m_dataFileName = dataFileName;
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
  if(!summaryData.checkVersion()) 
    {
      std::cout<<std::hex<<word[0]<<std::endl;
      std::cout<<word[1]<<std::endl;
      std::cout<<word[2]<<std::endl;
      std::cout<<word[3]<<std::dec<<std::endl;
      return false;
    }
  currentEvent->setEventSummaryData(summaryData);
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
      unsigned rawword = readRawWord();
      remainingBytes-=ANCILLARY_WORD_LENGTH;
      AncillaryWord OneWord;
      OneWord.setData(rawword);
      //      AncillaryWord OneWord = getNextDataWord();
      unsigned Header = OneWord.getHeader();
      QdcHeaderWord qdcHW;
      FadcHeaderWord fadcHW;
      ScalerHeaderWord scalerHW;
      
      //      std::cout<<"remainingBytes "<<remainingBytes<<std::endl;

      switch (Header)
	{
	case ANCILLARY_QDC_ID:
	  qdcHW.setData(rawword);
	  //	  std::cout<<"Check QDC Header: "<<qdcHW.checkHeader()<<" "<<qdcHW.getQcdFifo()<<std::endl;

	  for(int i=0; i < qdcHW.getQcdFifo(); i++)
	    {
	      QdcDataWord qdcDW;
	      qdcDW.setData(readRawWord());
	      currentEvent->appendQdcDataWord(qdcDW);
	    }
	  remainingBytes-=qdcHW.getQcdFifo()*ANCILLARY_WORD_LENGTH;
	  break;
	case ANCILLARY_FADC_ID:
	  fadcHW.setData(rawword);
	  //	  std::cout<<"Check FADC Header: "<<fadcHW.checkHeader()<<std::endl;
	  for(int i=0; i < fadcHW.getFadcFifo(); i++)
	    {
	      FadcDataWord fadcDW;
	      fadcDW.setData(readRawWord());
	      currentEvent->appendFadcDataWord(fadcDW);
	    }
	  remainingBytes-=fadcHW.getFadcFifo()*ANCILLARY_WORD_LENGTH;
	  break;
	case ANCILLARY_SCALER_ID:
	  scalerHW.setData(rawword);
	  //	  std::cout<<"Check SCALER Header: "<<scalerHW.checkHeader()<<" "<<scalerHW.getScalerFifo()<<std::endl;
	  for(int i=0; i < scalerHW.getScalerFifo(); i++)
	    {
	      ScalerDataWord scalerDW;
	      scalerDW.setData(readRawWord());
	      currentEvent->appendScalerDataWord(scalerDW);
	    }
	  remainingBytes-= scalerHW.getScalerFifo()*ANCILLARY_WORD_LENGTH;
	  break;
	default:
	  break;
	}
    }
  return true;
}
  
AdfEvent *AncillaryDataServer::getNextEvent()
{
  // Get the event header:
  AdfEvent *currentEvent = new AdfEvent();
  bool suxess;
  suxess = readEventHeader(currentEvent);
  if(!suxess) 
    {
      std::cout<<"failed to read the hevent header: "<<std::endl;
      currentEvent->getEventSummaryData().printEventHeader();
      return 0;
    }
  suxess = readEventData(currentEvent);
  if(!suxess) return 0;
  return currentEvent;
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

