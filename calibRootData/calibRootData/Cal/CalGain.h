#ifndef calibRootData_CalGain_h
#define calibRootData_CalGain_h

#include "TObject.h"
#include "TString.h"
#include "commonRootData/idents/calXtalId.h"
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
  class CalGain : public TObject  {
    
  public:

    /// Default constructor sets range, face to disallowed values
    CalGain(UInt_t packedId=0, Float_t val=0.0, Float_t sig=0.0) :
      m_id(packedId), m_gain(val), m_sig(sig) {}

    virtual ~CalGain() { };

  private:
    commonRootData::CalXtalId  m_id;
    Float_t  m_gain;
    Float_t  m_sig;

    ClassDef(calibRootData::CalGain, 1)
  };
}
#endif
