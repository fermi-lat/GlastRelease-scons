// $Header$

#include <string>
#include <vector>
#include "XmlBaseCnv.h"
// #include "XmlSAABoundaryCnv.h"
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

#include "CalibData/Nas/CalibSAABoundary.h"
#include "CalibData/CalibTime.h"
#include "xmlBase/Dom.h"

// Temporary.  Hope to find a better way to do this
#include "CalibData/CalibModel.h"

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

//static CnvFactory<XmlSAABoundaryCnv> s_factory;
//const  ICnvFactory& XmlSAABoundaryCnvFactory = s_factory;
DECLARE_CONVERTER_FACTORY(XmlSAABoundaryCnv);


XmlSAABoundaryCnv::XmlSAABoundaryCnv( ISvcLocator* svc) :
  XmlBaseCnv(svc, CLID_Calib_NAS_SAABoundary) { 
}


const CLID& XmlSAABoundaryCnv::objType() const {
  return CLID_Calib_NAS_SAABoundary;
}

const CLID& XmlSAABoundaryCnv::classID() {
  return CLID_Calib_NAS_SAABoundary;
}

 
// Create our specific object
StatusCode XmlSAABoundaryCnv::i_createObj(const DOMElement* element, 
                                    DataObject*& refpObject)
{
  using xmlBase::Dom;

  // Fetch quantities we need: the edges of the polygon
  DOMElement* polygon = Dom::findFirstChildByName(element, "polygon");
  if (polygon == 0) return StatusCode::FAILURE;
  DOMElement* vertex = Dom::findFirstChildByName(polygon, "vertex");
  if (vertex == 0) return StatusCode::FAILURE;

  while(vertex!=0){
      try {
	m_lat.push_back(Dom::getDoubleAttribute(vertex, "lat"));
	m_lon.push_back(Dom::getDoubleAttribute(vertex, "lon"));
      }  
      catch (xmlBase::DomException ex) {
	std::cerr << "From CalibSvc::XmlTest1Cnv::i_crateObj " << std::endl;
	std::cerr << ex.getMsg() << std::endl;
      }
      vertex=Dom::getSiblingElement(vertex);  
  };
    
  refpObject = new  
    CalibData::CalibSAABoundary(m_lat, m_lon ,*m_vstart, *m_vend,m_serNo);

  return StatusCode::SUCCESS;
}
