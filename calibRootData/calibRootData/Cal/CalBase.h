#ifndef calibRootData_CalBase_h
#define calibRootData_CalBase_h

#include "calibRootData/Base.h"
#include "calibRootData/Cal/CalDimension.h"

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
      Base(i, t), m_dim(nRow, nCol, nLayer, nXtal, nFace, nRange), 
      m_exact(false) { }
    virtual ~CalBase() { };

    const CalDimension* getCalDimension() const {return &m_dim;}
    bool getExact() const {return m_exact;}
  protected:

    CalDimension m_dim;
    bool_t m_exact;                       // true iff dimension info is exact
    ClassDef(calibRootData::CalBase, 1)
  };
}
#endif
