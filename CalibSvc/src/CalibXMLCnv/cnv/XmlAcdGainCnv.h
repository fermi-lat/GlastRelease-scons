// $Header$
#ifndef CalibSvc_XmlAcdGainCnv_h
#define CalibSvc_XmlAcdGainCnv_h

/** @class XmlAcdGainCnv

  Converter from xml to TCDS ACD gains class

  @author J. Bogart
*/
#include "XmlAcdBaseCnv.h"

template <class TYPE> class CnvFactory;

class XmlAcdGainCnv : public XmlAcdBaseCnv {

  /// Friend needed for instantiation
  friend class CnvFactory<XmlAcdGainCnv>;
public:
  const CLID& objType() const;
  static const CLID& classID();
protected:

  XmlAcdGainCnv(ISvcLocator* svcs);

  virtual ~XmlAcdGainCnv() {}       // most likely nothing to do 

  virtual StatusCode i_createObj(const DOMElement* element,
                                 DataObject*& refpObject);

};


#endif
