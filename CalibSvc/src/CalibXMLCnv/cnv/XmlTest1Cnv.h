// $Header$
#ifndef CalibData_XmlTest1Cnv_h
#define CalibData_XmlTest1Cnv_h

/** @class XmlTest1  

  Converter from xml to very simple calibration data-like class to be 
  used for testing calibration infrastructure

  @author J. Bogart
*/
#include "XmlBaseCnv.h"
#include <dom/DOM_Element.hpp>

class XmlTest1Cnv : public XmlBaseCnv {

  /// Friend needed for instantiation
  friend class CnvFactory<XmlTest1Cnv>;
public:
  const CLID& objType() const;
  static const CLID& classID();
protected:

  XmlTest1Cnv(ISvcLocator* svcs);

  virtual ~XmlTest1Cnv() {}       // most likely nothing to do 

  virtual StatusCode i_createObj(const DOM_Element& element,
                                 DataObject*& refpObject);

  
};


#endif
