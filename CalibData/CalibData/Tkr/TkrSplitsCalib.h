// $Header$
#ifndef CalibData_TkrSplitsCalib_h
#define CalibData_TkrSplitsCalib_h

#include "CalibData/Tkr/TkrBase.h"
#include "CalibData/CalibModel.h"
#include "CalibData/RangeBase.h"
#include <vector>

/** 
       @file CalibData/TkrSplits.h

       Defines structure of splits data: 2 numbers per Si uniplane to
       indicate which side reads out which strips

       @author  J. Bogart
*/

namespace CalibData {

  class TkrSplit : public RangeBase {
  public:
    TkrSplit(unsigned short low=0, unsigned short high=0) :
    m_low(low), m_high(high) {};
    ~TkrSplit() {};
    unsigned getLow() const {return m_low;}
    unsigned getHigh() const {return m_high;}

    // reimplemented from RangeBase
    virtual void update(RangeBase* other);

  private:
    unsigned short m_low;
    unsigned short m_high;
  };

  class TkrSplitsCalib : public TkrBase {
    friend class XmlTkrSplitsCnv;

  public:
    TkrSplitsCalib(unsigned nTowerRow, unsigned nTowerCol, 
                   unsigned nTray) : 
      TkrBase(nTowerRow, nTowerCol, nTray, 0, false) 
    {
      //      m_splits.reserve(m_finder->getSize());
      m_splits.resize(m_finder->getSize());
    }

    ~TkrSplitsCalib() {m_splits.clear();};


    // Reimplemented from DataObject
    virtual const CLID& clID() const {return classID(); }
    static const CLID&  classID();

    // Reimplemented from TkrBase
    virtual RangeBase* getChannel(const idents::TkrId& id, unsigned feChip=0);
    virtual RangeBase* getChannel(unsigned towerRow, unsigned towerCol,
                                  unsigned tray, bool top, unsigned feChip=0);
    
    virtual bool putChannel(RangeBase* data, const idents::TkrId& id, 
                            unsigned feChip=0);

    virtual bool putChannel(RangeBase* data, unsigned towerRow, 
                            unsigned towerCol, unsigned tray, 
                            bool top, unsigned feChip=0);

    virtual StatusCode update(CalibBase& other, MsgStream* log);
    
  private:
    std::vector<TkrSplit> m_splits;

  };

}
#endif
