// $Header$

#include <string>
#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/IOpaqueAddress.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IAddressCreator.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IConversionSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/GenericAddress.h"

#include "CalibSvc/ICalibXmlSvc.h"
#include "CalibSvc/ICalibMetaCnvSvc.h"

#include "CalibData/Nas/CalibLATAlignment.h"
#include "CalibData/CalibTime.h"
#include "xmlBase/Dom.h"

// Temporary.  Hope to find a better way to do this
#include "CalibData/CalibModel.h"
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

static CnvFactory<XmlLATAlignmentCnv> s_factory;
const  ICnvFactory& XmlLATAlignmentCnvFactory = s_factory;


XmlLATAlignmentCnv::XmlLATAlignmentCnv( ISvcLocator* svc) :
  XmlBaseCnv(svc, CLID_Calib_NAS_LATAlignment) { 
}


const CLID& XmlLATAlignmentCnv::objType() const {
  return CLID_Calib_NAS_LATAlignment;
}

const CLID& XmlLATAlignmentCnv::classID() {
  return CLID_Calib_NAS_LATAlignment;
}

 
// Create our specific object
StatusCode XmlLATAlignmentCnv::i_createObj(const DOMElement* element, 
                                           DataObject*& refpObject)
{
   using xmlBase::Dom;

  // Fetch quantities we need: the edges of the polygon
   DOMElement* angle = Dom::findFirstChildByName(element, "angle");
   if (angle == 0) return StatusCode::FAILURE;

   try {
     m_roll   = Dom::getDoubleAttribute(angle, "roll");
     m_pitch = Dom::getDoubleAttribute(angle, "pitch");
     m_yaw   = Dom::getDoubleAttribute(angle, "yaw");
   } catch (xmlBase::DomException ex) {
	std::cerr << "From CalibSvc::XmlLATAlignmentCnv::i_crateObj " << std::endl;
	std::cerr << ex.getMsg() << std::endl;
   }
    
   refpObject = new  
     CalibData::CalibLATAlignment(m_roll, m_pitch, m_yaw, *m_vstart, *m_vend,m_serNo);

   return StatusCode::SUCCESS;
}
