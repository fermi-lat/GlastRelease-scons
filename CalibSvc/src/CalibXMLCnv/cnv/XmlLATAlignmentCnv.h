// $Header$
#ifndef CalibSvc_XmlLATAlignment_h
#define CalibSvc_XmlLATAlignment_h

/** @class XmlTest1  

  Converter from xml to to a SAA calibration data object

  @author M. Ackermann
*/
#include "XmlBaseCnv.h"

template <class TYPE> class CnvFactory;

class XmlLATAlignmentCnv : public XmlBaseCnv {

  /// Friend needed for instantiation
  friend class CnvFactory<XmlLATAlignmentCnv>;
public:
  const CLID& objType() const;
  static const CLID& classID();
protected:

  XmlLATAlignmentCnv(ISvcLocator* svcs);

  virtual ~XmlLATAlignmentCnv() {}       // most likely nothing to do 

  virtual StatusCode i_createObj(const DOMElement* element,
                                 DataObject*& refpObject);
				 
  double m_roll;
  double m_pitch;
  double m_yaw;				 
   
};


#endif
