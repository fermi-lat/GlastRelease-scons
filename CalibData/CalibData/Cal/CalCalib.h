#include "CalibData/Cal/CalFinder.h"
#include "CalibData/CalibBase.h"

template <class CalCalibRange> class CalCalib : CalibData::CalibBase {
public: 

  inline CalCalib(const CalibData::CalFinder& finder); 
  inline ~CalCalib();

  inline CalCalib* getData(CalXtalId id, unsigned range, unsigned face=0) {
    unsigned ix = m_finder.findIx(id, range, face);
    if (ix < m_finder.getSize() ) { 
      return m_data[ix];
    }
    else return 0;
  }

  inline bool putRangeData(CalCalibRange *pData, CalXtalId id, 
                           unsigned range, unsigned face=0);

  // Following two functions are re-implemented from DataObject
  inline virtual const CLID& clID() const { classID(); }
    
  // We insist that the template parameter class have a static
  // function classID()
  static inline const CLID& classID() { return CalCalibRange::classID(); }


    // Re-implemented from CalibBase
  inline virtual void update(CalibData::CalibBase& other);

private: 
  CalCalibRange** m_data;
  const CalibData::CalFinder& m_finder;


};

inline CalCalib<CalCalibRange>::CalCalib(const CalibData::CalFinder& finder) : 
  m_finder(finder), m_data(0)                    {
  typedef CalCalibRange* CalCalibRange_ptr;
  
  unsigned s = m_finder.getSize();
  unsigned i;

  m_data = new CalCalibRange_ptr[s];
  for (i = 0; i < s; i++) m_data[i] = 0;
}

inline bool CalCalib<CalCalibRange>::putRangeData(CalCalibRange *pData, 
                                                  CalXtalId id, 
                                                  unsigned range, 
                                                  unsigned face=0) {

  unsigned ix = m_finder.findIx(id, range, face);
  if (ix < m_finder.getSize() ) {

    // CalCalibRange has to have a sensible copy constructor
    if (m_data[ix]) delete m_data[ix];
    m_data[ix] = new CalCalibRange(*pData);
    return true;
  }
  return false;
}
      

inline CalCalib<CalCalibRange>::~CalCalib() {
  unsigned i;
  unsigned s = m_finder.getSize();
  for (i = 0; i < s; i++) {
    delete m_data[i];
  }
  delete [] m_data;
}

inline virtual void CalCalib<CalCalibRange>::update(CalibData::CalibBase& 
                                                    other) {
  CalCalib& other1 = dynamic_cast<CalCalib<CalCalibRange>& >(other);
  CalibBase::update(other1);

  // Ought to be the case that finder is already OK; just need
  // to update data
  unsigned i;
  unsigned s = m_finder.getSize();
  for (i = 0; i < s; i++) {
    if (m_data[i]) delete m_data[i];
    m_data[i] = new CalCalibRange(*other.m_data[i]);
  }
}
  


