#ifndef calibRootData_CalBase_h
#define calibRootData_CalBase_h

#include "../Base.h"
#include "../Cal/CalDimension.h"

// Don't need these since they're already in Base
// #include "TObject.h"
// #include "TString.h"
/*  Might need some more ROOT includes
#include "TRef.h"
#include "TRefArray.h"
*/

/** @class calibRootData::CalBase
* @brief Base class for (almost) all Cal calibration ROOT classes.  
* 
* Those representing hot/dead channels probably won't use this base.
*
* @author Joanne Bogart
*  
* $Header$
*/
namespace calibRootData {
  class CalBase : public Base  {
    
  public:


    CalBase(TString i = TString(), TString t = TString(),
            unsigned nRow=0, unsigned nCol=0, unsigned nLayer=0, 
            unsigned nXtal=0, unsigned nFace=0, unsigned nRange=0) : 
      Base(i, t), m_dim(nRow, nCol, nLayer, nXtal, nFace, nRange, false)
      { }
    virtual ~CalBase() { };

    const CalDimension* getCalDimension() const {return &m_dim;}
    bool getExact() const {return getCalDimension()->isExact();}

    void initialize(unsigned nRow, unsigned nCol, unsigned nLayer, 
                    unsigned nXtal, unsigned nFace, unsigned nRange,
                    unsigned nDacCol) 
    {
      m_dim.initialize(nRow, nCol, nLayer, nXtal, nFace, nRange, nDacCol);
    }
    void setExact(bool exact) {m_dim.m_exact = exact;}
  protected:

    CalDimension m_dim;
    ClassDef(calibRootData::CalBase, 1)
  };
}
#endif
