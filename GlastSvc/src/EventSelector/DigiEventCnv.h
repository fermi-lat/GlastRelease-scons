// DigiEventCnv.h: interface for the DigiEventCnv class.
//
//////////////////////////////////////////////////////////////////////

#ifndef CNV_DIGIEVENTCNV_H 
#define CNV_DIGIEVENTCNV_H 1

// Include files
#include "ItemCnv.h"

// Forward declarations
class IDataDirectory;
class ISvcLocator;
class DigiEvent;
// Abstract factory to create the converter
template <class T> class CnvFactory;

//!Used to make Digi branch of the TDS work.
/*! Right now this converter does nothing particularly useful other than
    make the Digi branch of the TDS work.
*/
class DigiEventCnv : public ItemCnv<DigiEvent>      {
  friend class CnvFactory<DigiEventCnv>;
protected:
  /// Standard Constructor
  DigiEventCnv(ISvcLocator* svc);
  /// Standard Destructor
  virtual ~DigiEventCnv();
  /// Update object from somewhere
  virtual StatusCode updateObj(int* data, DigiEvent* pObject);
};
#endif // CNV_DIGIEVENTCNV_H
