// $Header$
#ifndef CalibData_TkrSplitsCalib_h
#define CalibData_TkrSplitsCalib_h

#include "CalibData/Tkr/OldTkrBase.h"
#include "CalibData/CalibModel.h"
#include "CalibData/RangeBase.h"
#include <vector>

/** 
       @file CalibData/TkrSplits.h

       Defines structure of splits data: 2 numbers per Si uniplane to
       indicate which side reads out which strips

       @author  J. Bogart
*/

class XmlTkrSplitsCnv;

namespace CalibData {

  class TkrSplit : public RangeBase {
  public:
    TkrSplit(unsigned short low=0, unsigned short high=0) :
    m_low(low), m_high(high) {};
    virtual ~TkrSplit() {};
    unsigned getLow() const {return m_low;}
    unsigned getHigh() const {return m_high;}

    // reimplemented from RangeBase
    virtual void update(RangeBase* other);

  private:
    unsigned short m_low;
    unsigned short m_high;
  };

  class TkrSplitsCalib : public OldTkrBase {
    friend class ::XmlTkrSplitsCnv;

  public:
    TkrSplitsCalib(unsigned nTowerRow, unsigned nTowerCol, 
                   unsigned nTray) : 
      OldTkrBase(nTowerRow, nTowerCol, nTray, false) 
    {
      //      m_splits.reserve(m_finder->getSize());
      m_splits.resize(m_finder->getSize());
    }

    ~TkrSplitsCalib() {m_splits.clear();};


    // Reimplemented from DataObject
    virtual const CLID& clID() const {return classID(); }
    static const CLID&  classID();

    // Reimplemented from OldTkrBase
    virtual RangeBase* getChannel(const idents::TkrId& id);
    virtual RangeBase* getChannel(unsigned towerRow, unsigned towerCol,
                                  unsigned tray, bool top);
    
    virtual bool putChannel(RangeBase* data, const idents::TkrId& id);

    virtual bool putChannel(RangeBase* data, unsigned towerRow, 
                            unsigned towerCol, unsigned tray, 
                            bool top);

    virtual StatusCode update(CalibBase& other, MsgStream* log);
    
  private:
    std::vector<TkrSplit> m_splits;

  };

}
#endif
