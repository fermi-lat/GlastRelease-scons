#ifndef calibRootData_CalSerial_h
#define calibRootData_CalSerial_h

#include "TObject.h"
#include "TString.h"
/*  Might need some more ROOT includes
#include "TRef.h"
#include "TRefArray.h"
*/

/** @class calibRootData::CalSerial
* @brief Keep track of cal serial numbers versus grid bay
* 
* Note grid bay numbering, as viewed from above (TKR side), looks like this:
*
*       12   13   14   15
*        8    9   10   11
*        4    5    6    7
*        0    1    2    3
*
* @author Joanne Bogart
*  
* $Header$
*/
namespace calibRootData {
  class CalSerial : public TObject {
  public:
    
    CalSerial() {}

    ~CalSerial() { };

    class GridSerial {
    public:
      UInt_t gridBay;
      TString serial;
    };
    void addEntry(calibRootData::CalSerial::GridSerial entry) {m_entries.push_back(entry);}
    const std::vector<calibRootData::CalSerial::GridSerial>& getEntries() const {return m_entries;}
  private:

      std::vector<calibRootData::CalSerial::GridSerial> m_entries;
  
    ClassDef(calibRootData::CalSerial, 1)
  };
}
#endif
