// $Header$
#ifndef CalibSvc_XmlCalLightAttCnv_h
#define CalibSvc_XmlCalLightAttCnv_h

/** @class XmlCalLightAttCnv

  Converter from xml to TCDS CAL light attenuation class

  @author J. Bogart
*/
#include "XmlCalBaseCnv.h"
#include <xercesc/dom/DOM_Element.hpp>

template <class TYPE> class CnvFactory;

class XmlCalLightAttCnv : public XmlCalBaseCnv {

  /// Friend needed for instantiation
  friend class CnvFactory<XmlCalLightAttCnv>;
public:
  const CLID& objType() const;
  static const CLID& classID();
protected:

  XmlCalLightAttCnv(ISvcLocator* svcs);

  virtual ~XmlCalLightAttCnv() {}       // most likely nothing to do 

  virtual StatusCode i_createObj(const DOM_Element& element,
                                 DataObject*& refpObject);

};


#endif
