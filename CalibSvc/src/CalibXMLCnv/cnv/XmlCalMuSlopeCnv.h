// $Header$
#ifndef CalibSvc_XmlCalMuSlopeCnv_h
#define CalibSvc_XmlCalMuSlopeCnv_h

/** @class XmlCalMuSlopeCnv

  Converter from xml to TCDS CAL mu slope class

  @author J. Bogart
*/
#include "XmlBaseCnv.h"
#include <xercesc/dom/DOM_Element.hpp>

template <class TYPE> class CnvFactory;

class XmlCalMuSlopeCnv : public XmlBaseCnv {

  /// Friend needed for instantiation
  friend class CnvFactory<XmlCalMuSlopeCnv>;
public:
  const CLID& objType() const;
  static const CLID& classID();
protected:

  XmlCalMuSlopeCnv(ISvcLocator* svcs);

  virtual ~XmlCalMuSlopeCnv() {}       // most likely nothing to do 

  virtual StatusCode i_createObj(const DOM_Element& element,
                                 DataObject*& refpObject);

};


#endif
