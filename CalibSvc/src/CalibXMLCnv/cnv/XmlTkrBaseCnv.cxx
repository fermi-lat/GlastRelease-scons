// $Header$
#include "GaudiKernel/MsgStream.h"
#include "XmlTkrBaseCnv.h"
#include "xml/Dom.h"
#include "idents/TkrId.h"


XmlTkrBaseCnv::XmlTkrBaseCnv(ISvcLocator* svc, const CLID& clid) :
  XmlBaseCnv(svc, clid), m_nChannel(0) {}

StatusCode XmlTkrBaseCnv::readDimension(const DOM_Element& docElt, 
                                        unsigned& nRow, unsigned& nCol, 
                                        unsigned& nTray, unsigned& nChip)
{
  using xml::Dom;

  MsgStream log(msgSvc(), "XmlTkrBaseCnv" );
  DOM_Element dimElt = Dom::findFirstChildByName(docElt, "dimension");
  if (dimElt == DOM_Element()) return StatusCode::FAILURE;

  try {
    nRow = Dom::getIntAttribute(dimElt, "nBayRow");
    nCol = Dom::getIntAttribute(dimElt, "nBayCol");
    nTray = Dom::getIntAttribute(dimElt, "nTray");
    nChip = 0;
    std::string chipStr = Dom::getAttribute(dimElt, "nChip");
    if (chipStr.size() > 0) {
      nChip = Dom::getIntAttribute(dimElt, "nChip");
    }
  }
  catch (xml::DomException ex) {
    std::cerr << "From CalibSvc::XmlTkrBaseCnv::readDimension" << std::endl;
    std::cerr << ex.getMsg() << std::endl;
    throw ex;
  }

  return StatusCode::SUCCESS;
}

