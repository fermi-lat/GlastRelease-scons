#ifndef calibRootData_CalGain_h
#define calibRootData_CalGain_h

#include "TObject.h"
#include "commonRootData/idents/CalXtalId.h"

/** @class calibRootData::CalGain
* @brief  Gain information for a single channel (includes channel identifier)
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

    const CalXtalId& getId() const {return m_id;}
    const float getGain() const {return m_gain;}
    const float getSig() const {return m_sig;}
    

  private:
    CalXtalId  m_id;
    Float_t  m_gain;
    Float_t  m_sig;

    ClassDef(calibRootData::CalGain, 1)
  };
}
#endif
