// $Header$

/// @file AcdCalibBase
/// @author J. Bogart
#ifndef CalibData_AcdCalibBase_h
#define CalibData_AcdCalibBase_h

#include <vector>
#include "CalibData/CalibBase.h"
#include "idents/AcdId.h"

namespace CalibData {
  class AcdFinder;
  class RangeBase;

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
                 unsigned nNA=11, unsigned nPmt=2);

    virtual ~AcdCalibBase();

    /** 
        Pick out calibration data associated with a particular tile, pmt
     */
    virtual RangeBase* getPmt(idents::AcdId id, unsigned pmt=0);

    bool putPmt(idents::AcdId id, unsigned pmt, RangeBase* data);

    /* Do we need this version?  I hope not - doesn't cover NAs */
    bool putPmt(unsigned face, unsigned row, unsigned col, 
                  unsigned pmt, RangeBase* data);


    virtual const CLID& clID() const = 0;     // must be overridden
    static const CLID& classID();   // shouldn't get called

    // Maybe won't need to be virtual after all
    virtual StatusCode update(CalibBase& other, MsgStream* log);
    
  protected:
    AcdFinder* m_finder;

    std::vector<RangeBase* > m_pmts;  // tiles and ribbons
    std::vector<RangeBase* > m_NAs;   // no detector

  private:
    static const CLID noCLID;

    /** Due to bug in gcc, gdb can't find symbols in constructors.  This
        method is called by the constructor and does most of the work
    */
    void cGuts(unsigned nFace, unsigned nRow, unsigned nCol, 
               unsigned nNA, unsigned nPmt);

  };

}  
#endif
