// $Header$
#ifndef CalibSvc_XmlCalPedCnv_h
#define CalibSvc_XmlCalPedCnv_h

/** @class XmlCalPed1  

  Converter from xml to TCDS CAL pedestals class

  @author J. Bogart
*/
#include "XmlBaseCnv.h"
#include <xercesc/dom/DOM_Element.hpp>

template <class TYPE> class CnvFactory;

class XmlCalPedCnv : public XmlBaseCnv {

  /// Friend needed for instantiation
  friend class CnvFactory<XmlCalPedCnv>;
public:
  const CLID& objType() const;
  static const CLID& classID();
protected:

  XmlCalPedCnv(ISvcLocator* svcs);

  virtual ~XmlCalPedCnv() {}       // most likely nothing to do 

  virtual StatusCode i_createObj(const DOM_Element& element,
                                 DataObject*& refpObject);

};


#endif
