#ifndef calibRootData_Base_h
#define calibRootData_Base_h

#include "TObject.h"
#include "TString.h"
/*  Might need some more ROOT includes
#include "TRef.h"
#include "TRefArray.h"
*/

/** @class BaseCalibRoot
* @brief Base class for all calibration ROOT classes.
*
* @author Joanne Bogart
*  
* $Header$
*/
namespace calibRootData {
  class Base : public TObject {
    
  public:

    /// return true if arg. matches internal instrument name
    bool isInstrument(const TString instrument) const {
      return (instrument == m_instrument);
    }

    bool isType(const TString calibType) const {
      return (calibType == m_calibType); 
    }

    // Should these arguments default to TString equivalent of empty?
    Base(TString i, TString t) : m_instrument(i), m_calibType(t) { };

    void setInstrument(TString i) {
      m_instrument = i; 
    }
    virtual ~Base() { };

  protected:
    TString m_instrument;
    TString m_calibType;

    // Do we want anything else here?  
    // Something else of general interest is format version, but ROOT has
    // its own way to deal with schema evolution.  
    //
    // Other stuff which goes in base TDS class for calib comes from metadata.
    // (instrument and calibType can also be found there, and in fact those 
    // are the crucial ones, not the ones internal to file.  Probably should
    // check to see if they match, and complain if they don't).

    ClassDef(calibRootData::Base, 1)

  };
}
#endif
