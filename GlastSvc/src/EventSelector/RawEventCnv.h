// RawEventCnv.h: interface for the RawEventCnv class.
//
//////////////////////////////////////////////////////////////////////

#ifndef CNV_RAWEVENTCNV_H 
#define CNV_RAWEVENTCNV_H 1

// Include files
#include "ItemCnv.h"

// Forward declarations
class IDataDirectory;
class ISvcLocator;
class RawEvent;
// Abstract factory to create the converter
template <class T> class CnvFactory;

//!Used to make Raw branch of the TDS work.
/*! Right now this converter does nothing particularly useful other than
    make the Raw branch of the TDS work.
*/
class RawEventCnv : public ItemCnv<RawEvent>      {
  friend class CnvFactory<RawEventCnv>;
protected:
  /// Standard Constructor
  RawEventCnv(ISvcLocator* svc);
  /// Standard Destructor
  virtual ~RawEventCnv();
  /// Update object from somewhere
  virtual StatusCode updateObj(int* data, RawEvent* pObject);
};
#endif // CNV_RAWEVENTCNV_H
