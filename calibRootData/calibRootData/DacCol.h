#ifndef CALIBROOTDATA_DACCOL_H
#define CALIBROOTDATA_DACCOL_H

#include "TObject.h"
#include <vector>
/*  Might need some more ROOT includes
#include "TRef.h"
#include "TRefArray.h"
*/

/** @class calibRootData::DacCol
* @brief Utility calibration class to describe a collection of dac settings
*
* @author Joanne Bogart
*  
* $Header$
*/
namespace calibRootData {
  class DacCol : public TObject {
    
  public:

    // Should these arguments default to TString equivalent of empty?
    DacCol() {m_dacs.clear();};

    // Also need copy constructor

    virtual ~DacCol() {m_dacs.clear(); };

    std::vector<UShort_t>& getDacs() const {return m_dacs;}

    void append(UShort_t value) {m_dacs.push_back(value);}

  private:
    
    std::vector<UShort_t> m_dacs;

    // Do we want anything else here?  
    // Something else of general interest is format version, but ROOT has
    // its own way to deal with schema evolution.  
    //
    // Other stuff which goes in base TDS class for calib comes from metadata.
    // (instrument and calibType can also be found there, and in fact those 
    // are the crucial ones, not the ones internal to file.  Probably should
    // check to see if they match, and complain if they don't).

    ClassDef(calibRootData::DacCol, 1)

  };
}
#endif
