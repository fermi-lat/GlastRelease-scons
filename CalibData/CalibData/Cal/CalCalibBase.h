// $Header$

/// @file CalCalibBase
/// @author J. Bogart
#ifndef CalibData_CalCalibBase_h
#define CalibData_CalCalibBase_h

#include <vector>
#include "CalibData/CalibBase.h"
#include "idents/CalXtalId.h"

namespace CalibData {
  class CalFinder;
  class RangeBase;

  /**
       Base class for calorimeter calibration data, at least for those
       types with a fixed amount of data per crystal-range-face

       This class keeps a (pointer to a) vector of pointers to the
       individual per-range datasets and a pointer to a helper class,
       CalFinder, which knows how to compute indices.
  */
  class CalCalibBase : public CalibBase {

  public:
    CalCalibBase(unsigned nTowerRow=4, unsigned nTowerCol=4, unsigned nLayer=8, 
                 unsigned nXtal=12, unsigned nFace=2, unsigned nRange=4);
    virtual ~CalCalibBase();

    RangeBase* getRange(idents::CalXtalId id, 
                        unsigned range=0, unsigned face=0);

    bool putRange(idents::CalXtalId id, unsigned range, unsigned face, 
                  RangeBase* data);

    bool putRange(unsigned towerRow, unsigned towerCol, 
                  unsigned layer, unsigned xtal, unsigned range,
                  unsigned face, RangeBase* data);

    virtual const CLID& clID() const = 0;     // must be overridden
    static const CLID& classID();   // shouldn't get called

    // Maybe won't need to be virtual after all
    virtual StatusCode update(CalibBase& other, MsgStream* log);
    
  protected:
    CalFinder* m_finder;
    //    std::vector<RangeBase* >* m_pR;
    std::vector<RangeBase* > m_ranges;
  private:
    static const CLID noCLID;
  };

}  
#endif
