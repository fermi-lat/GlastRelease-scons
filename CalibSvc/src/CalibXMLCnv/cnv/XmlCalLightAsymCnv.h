// $Header$
#ifndef CalibSvc_XmlCalLightAsymCnv_h
#define CalibSvc_XmlCalLightAsymCnv_h

/** @class XmlCalLightAsymCnv

  Converter from xml to TCDS CAL LightAtt class

  @author J. Bogart
*/
#include "XmlCalBaseCnv.h"

template <class TYPE> class CnvFactory;

class XmlCalLightAsymCnv : public XmlCalBaseCnv {

  /// Friend needed for instantiation
  friend class CnvFactory<XmlCalLightAsymCnv>;
public:
  const CLID& objType() const;
  static const CLID& classID();
protected:

  XmlCalLightAsymCnv(ISvcLocator* svcs);

  virtual ~XmlCalLightAsymCnv() {}       // most likely nothing to do 

  virtual StatusCode i_createObj(const DOMElement* element,
                                 DataObject*& refpObject);

};


#endif
