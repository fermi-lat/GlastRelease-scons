// $Header$
#ifndef RootCalGainCnv_h
#define RootCalGainCnv_h

/** @class RootCalGainCnv

  Converter from Root to TCDS CAL gains class

  @author J. Bogart
*/
#include "RootCalBaseCnv.h"


template <class TYPE> class CnvFactory;

class RootCalGainCnv : public RootCalBaseCnv {

  /// Friend needed for instantiation
  friend class CnvFactory<RootCalGainCnv>;
public:
  const CLID& objType() const;
  static const CLID& classID();

  // For writing
  virtual StatusCode createRoot(const std::string& fname, 
                                CalibData::CalibBase* pTDSObj);


protected:

  /**
     Given a pointer to a TDS object which can be cast to "our" type, fill
     in corresponding information in the corresponding root class

     @param pTDSObj   Pointer to tds object to be converted
     @param pRootObj  Pointer to destination root object

  */
  virtual StatusCode fillRoot(CalibData::CalibBase* pTDSObj, 
                              TObject* pRootObj);

  RootCalGainCnv(ISvcLocator* svcs);

  virtual ~RootCalGainCnv() {}       // most likely nothing to do 

  virtual StatusCode i_createObj (const std::string& fname,
                                  DataObject*& refpObject);
  /*  virtual StatusCode i_createObj(TObject* pRootObj,       
      DataObject*& refpObject);  */

};


#endif
