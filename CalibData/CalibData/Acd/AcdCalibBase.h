// $Header$

/// @file CalCalibBase
/// @author J. Bogart
#ifndef CalibData_CalCalibBase_h
#define CalibData_CalCalibBase_h

#include <vector>
#include "CalibData/CalibBase.h"
#include "idents/AcdId.h"

namespace CalibData {
  class AcdFinder;
  class RangeBase;
  //  class DacCol;

  /**
       Base class for ACD calibration data, at least for those
       types with a fixed amount of data per tile-pmt-range

       This class keeps a (pointer to a) vector of pointers to the
       individual per-range datasets and a pointer to a helper class,
       AcdFinder, which knows how to compute indices.
  */
  class AcdCalibBase : public CalibBase {

  public:
    AcdCalibBase(unsigned nFace=5, unsigned nRow=5, unsigned nCol=5, 
                 unsigned nPmt=2, unsigned nRange=2,
                 unsigned nDacCol=0);
    virtual ~AcdCalibBase();

    /** 
        Pick out calibration data associated with a particular crystal,
        face, range.
        May need to be overridden in case same data should be associated
        with more than one range (e.g., light asym)
     */
    virtual RangeBase* getRange(idents::AcdId id, 
                                unsigned pmt=0, unsigned range=0);

    bool putRange(idents::AcdId id, unsigned pmt, unsigned range, 
                  RangeBase* data);

    bool putRange(unsigned face, unsigned row, unsigned col, 
                  unsigned pmt, unsigned range, RangeBase* data);

    // bool putDacCol(unsigned range, DacCol* dacs);

    //    DacCol* getDacCol(unsigned range);

    virtual const CLID& clID() const = 0;     // must be overridden
    static const CLID& classID();   // shouldn't get called

    // Maybe won't need to be virtual after all
    virtual StatusCode update(CalibBase& other, MsgStream* log);
    
  protected:
    AcdFinder* m_finder;
    //    std::vector<RangeBase* >* m_pR;
    std::vector<RangeBase* > m_ranges;
    //    std::vector<DacCol*> m_dacCols; // often there aren't any
  private:
    static const CLID noCLID;

    /** Due to bug in gcc, gdb can't find symbols in constructors.  This
        method is called by the constructor and does most of the work
    */
    void cGuts(unsigned nFace, unsigned nRow, unsigned nCol, 
               unsigned nPmt, unsigned nRange /*, unsigned nDacCol */);

  };

}  
#endif
