// $Header$
#ifndef CalibSvc_XmlCalGainCnv_h
#define CalibSvc_XmlCalGainCnv_h

/** @class XmlCalGainCnv

  Converter from xml to TCDS CAL gains class

  @author J. Bogart
*/
#include "XmlCalBaseCnv.h"

template <class TYPE> class CnvFactory;

class XmlCalGainCnv : public XmlCalBaseCnv {

  /// Friend needed for instantiation
  friend class CnvFactory<XmlCalGainCnv>;
public:
  const CLID& objType() const;
  static const CLID& classID();
protected:

  XmlCalGainCnv(ISvcLocator* svcs);

  virtual ~XmlCalGainCnv() {}       // most likely nothing to do 

  virtual StatusCode i_createObj(const DOMElement* element,
                                 DataObject*& refpObject);

};


#endif
