// $Header$

/// @file CalCalibBase
/// @author J. Bogart
#ifndef CalibData_CalCalibBase_h
#define CalibData_CalCalibBase_h

#include <vector>
#include "CalibData/CalibBase.h"
#include "CalibData/Cal/CalFinder.h"
#include "idents/CalXtalId.h"

namespace CalibData {
  class RangeBase;
  class DacCol;
  class Xpos;

  /**
       Base class for calorimeter calibration data, at least for those
       types with a fixed amount of data per crystal-range-face

       This class keeps a (pointer to a) vector of pointers to the
       individual per-range datasets and a pointer to a helper class,
       CalFinder, which knows how to compute indices.
  */
  class CalCalibBase : public CalibBase {

  public:
    CalCalibBase(unsigned nTowerRow=4, unsigned nTowerCol=4, unsigned nLayer=8,                  unsigned nXtal=12, unsigned nFace=2, unsigned nRange=4,
                 unsigned nDacCol=0, unsigned nXpos = 0);
    virtual ~CalCalibBase();

    /** 
        Pick out calibration data associated with a particular crystal,
        face, range.
        May need to be overridden in case same data should be associated
        with more than one range (e.g., light asym)
     */
    virtual RangeBase* getRange(idents::CalXtalId id, 
                                unsigned range=0, unsigned face=0);

    virtual RangeBase* getRange(unsigned towerRow, unsigned towerCol,
                                unsigned layer, unsigned xtal,
                                unsigned range=0, unsigned face=0);

    // **NEW** 
    bool putRange(idents::CalXtalId id, RangeBase* data) {
      return putRange(id, 0, 0, data);
    }

    bool putRange(idents::CalXtalId id, 
                  unsigned range, unsigned face, RangeBase* data);

    bool putRange(unsigned towerRow, unsigned towerCol, 
                  unsigned layer, unsigned xtal, unsigned range,
                  unsigned face, RangeBase* data);

    bool putDacCol(unsigned range, DacCol* dacs);

    DacCol* getDacCol(unsigned range) const;

    bool putXpos(Xpos* pos);

    Xpos* getXpos() const {return m_xpos;}

    // Get dimensioning information; needed when transforming to 
    // permanent storage
    /// Get # tower rows
    unsigned getNTowerRow() const {return m_finder->getNTowerRow();}

    /// Get # tower columns
    unsigned getNTowerCol() const {return m_finder->getNTowerCol();}

    /// Get #  layers
    unsigned getNLayer() const {return m_finder->getNLayer();}

    /// Get # crystals/layer
    unsigned getNXtal() const {return m_finder->getNXtal();}

    /// Get # face relevant for this calibration
    unsigned getNFace() const {return m_finder->getNFace();}

    /// Get # ranges relevant for this calibration
    unsigned getNRange() const {return m_finder->getNRange();}

    /// Get # dac setting collections
    unsigned getNDacCol() const {return m_finder->getNDacCol();}

    // Get # xpos collections 
    unsigned getNXpos() const {return m_finder->getNXpos();}

    virtual const CLID& clID() const = 0;     // must be overridden
    static const CLID& classID();   // shouldn't get called

    // Maybe won't need to be virtual after all
    virtual StatusCode update(CalibBase& other, MsgStream* log);
    
  protected:
    CalFinder* m_finder;
    //    std::vector<RangeBase* >* m_pR;
    std::vector<RangeBase* > m_ranges;
    std::vector<DacCol*> m_dacCols;              // often there aren't any
    Xpos* m_xpos;                                // usually 0
  private:
    static const CLID noCLID;

    /** Due to bug in gcc, gdb can't find symbols in constructors.  This
        method is called by the constructor and does most of the work
    */
    void cGuts(unsigned nTowerRow, unsigned nTowerCol, 
               unsigned nLayer, unsigned nXtal, unsigned nFace, 
               unsigned nRange, unsigned nDacCol, unsigned nXpos);

  };

}  
#endif
