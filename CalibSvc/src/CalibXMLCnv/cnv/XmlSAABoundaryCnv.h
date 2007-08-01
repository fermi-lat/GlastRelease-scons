// $Header$
#ifndef CalibSvc_XmlSAABoundary_h
#define CalibSvc_XmlSAABoundary_h

/** @class XmlTest1  

  Converter from xml to to a SAA calibration data object

  @author M. Ackermann
*/
#include "XmlBaseCnv.h"
#include <vector>

template <class TYPE> class CnvFactory;

class XmlSAABoundaryCnv : public XmlBaseCnv {

  /// Friend needed for instantiation
  friend class CnvFactory<XmlSAABoundaryCnv>;
public:
  const CLID& objType() const;
  static const CLID& classID();
protected:

  XmlSAABoundaryCnv(ISvcLocator* svcs);

  virtual ~XmlSAABoundaryCnv() {}       // most likely nothing to do 

  virtual StatusCode i_createObj(const DOMElement* element,
                                 DataObject*& refpObject);
  
  std::vector<double> m_lat;
  std::vector<double> m_lon;
   
};


#endif
