#ifndef CNV_DIGIEVENTCNV_H 
#define CNV_DIGIEVENTCNV_H 1

#include "BaseCnv.h"

class IDataDirectory;
class ISvcLocator;
class DigiEvent;

extern const CLID& CLID_DigiEvent;

// Abstract factory to create the converter
template <class T> class CnvFactory;

/** @class DigiEventCnv
 * @brief Right now this converter does nothing particularly useful other than
 * make the Digi branch of the TDS work.
 *
 * $Header$
*/
class DigiEventCnv : public BaseCnv {
  friend class CnvFactory<DigiEventCnv>;

public: 
  static const CLID& classID()   
  {
    return CLID_DigiEvent;
  }

protected:
  DigiEventCnv(ISvcLocator* svc);
  virtual ~DigiEventCnv() { };
  virtual StatusCode updateObj(int* data, DigiEvent* pObject);
};
#endif // CNV_DIGIEVENTCNV_H
