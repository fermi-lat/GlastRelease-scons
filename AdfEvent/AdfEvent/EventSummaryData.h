#ifndef EVENTSUMMARYDATA_HH
#define EVENTSUMMARYDATA_HH

#include <iostream>

/**
* @class EventSummaryData
* @brief Base class describing the event header words from the ancillary detectors.

* Based on ancillary data version 3.01. The basic structure is as follows:
* 31   20  19                                0
* -------  -----------------------------------
* version        ANCILLARY_EVENT_HEADER
* 31                                         0
* --------------------------------------------
*            event_length (bytes)
* 31              28  27                     0
* ------------------  ------------------------
* ANCILLARY_EVENT_ID        event_number
* 31              28  27        12  11       0
* ------------------  ------------  ----------
* ANCILLARY_EVENT_ID  spill_number  spill_size
*/


namespace AncillaryData {
  const unsigned ANCILLARY_DATA_VERSION = 301;
  const unsigned ANCILLARY_EVENT_HEADER = 0xf1030;
  const unsigned ANCILLARY_EVENT_ID     = 12;
  class EventSummaryData {
  public:
    EventSummaryData()  {m_word = 0;} 
    ~EventSummaryData() {;}//if (m_word) delete[] m_word;} 
    void setData(const unsigned int *word) 
      {
	m_word = new unsigned int[4];
	unsigned int i;
	for (i=0;i<4;i++)
	  m_word[i] = word[i];
	
	m_version     = (word[0] >> 20) & 0xfff;
	m_eventLength =  word[1];
	m_eventNumber = word[2] & 0xfffffff;
	m_spillNumber = (word[3] >> 12) & 0xffff;
	m_spillSize   = word[3] & 0xfff;
      }	
    
    
    unsigned getVersion()        const {return m_version;}
    unsigned getEventLength()    const {return m_eventLength;}
    unsigned getEventNumber()    const {return m_eventNumber;}
    unsigned getSpillNumber()    const {return m_spillNumber;}
    unsigned getSpillSize()      const {return m_spillSize;}
    /*
      unsigned int getVersion()     const {return unsigned((*m_word[0] >> 20) & 0xfff);}
      unsigned int getEventLength() const {return unsigned(*m_word[1]);}
      unsigned int getEventNumber() const {return unsigned(*m_word[2] & 0xfffffff);}
      unsigned int getSpillNumber() const {return unsigned((*m_word[3] >> 12) & 0xffff);}
      unsigned int getSpillSize()   const {return unsigned(*m_word[3] & 0xfff);}
    */
    void printEventHeader()
      {
	std::cout << "#-----------------------------------------------" << std::endl;
	std::cout << "Version :             " << getVersion()           << std::endl;
	std::cout << "Event length (bytes) :" << getEventLength()       << std::endl;
	std::cout << "Event number :        " << getEventNumber()       << std::endl;
	std::cout << "Spill number :        " << getSpillNumber()       << std::endl;
	std::cout << "Spill size :          " << getSpillSize()         << std::endl;
	std::cout << "#-----------------------------------------------" << std::endl;
      }
    bool checkVersion()              const {return (getVersion() == ANCILLARY_DATA_VERSION);}
    /*
      bool checkHeader()               const {return ((m_word[0] & 0xfffff) == ANCILLARY_EVENT_HEADER) &
      (((m_word[2] >> 28) & 0xf) == ANCILLARY_EVENT_ID) &
      (((m_word[3] >> 28) & 0xf) == ANCILLARY_EVENT_ID);}
    */
  private:
    unsigned int *m_word;
    unsigned m_version;
    unsigned m_eventLength;
    unsigned m_eventNumber;
    unsigned m_spillNumber;
    unsigned m_spillSize;
  };
 
}
#endif

