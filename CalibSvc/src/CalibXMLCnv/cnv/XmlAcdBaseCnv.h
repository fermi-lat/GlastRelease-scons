// $Header$
#ifndef CalibData_XmlAcdBaseCnv_h
#define CalibData_XmlAcdBaseCnv_h

/** @class XmlAcdBaseCnv 

  Base class for ACD calibration converters from XML files to TCDS.
  All such converters need to do certain things, which are
  handled here.  Methods common to *all* calibrations are in the
  base class XmlBaseCnv

  @author J. Bogart
*/

#include "XmlBaseCnv.h"
#include "idents/AcdId.h"
class XmlAcdBaseCnv : public XmlBaseCnv {

public:

  XmlAcdBaseCnv(ISvcLocator* svc, const CLID& clid);

  virtual ~XmlAcdBaseCnv() {};


  /** Convenience routine used by most ACD calibration types, which
   have a <dimension> element describing how the remainder of the
   data is laid out.
   @a nDet  number of tiles + ribbons
   @a nFace actually tile faces + ribbon orientations (7 for LAT)
   @a nRow  max row value to appear (5 for LAT)
   @a nCol  max col value or ribbon number (5 for LAT)
   @a nPmt  always = 2
   @a nNA   max number of unconnected electronics channels 
  */
  StatusCode readAcdDimension(const DOMElement* docElt, unsigned& nDet,
                              unsigned& nFace, unsigned& nRow,
                              unsigned& nCol, unsigned& nNA,
                              unsigned& nPmt);


  /// Another one to find first ped, gain, etc. (depending on calib. type)
  DOMElement* findFirstPmt(const DOMElement* docElt);

  /// Still another one to navigate XML file and find next set of pmt data
  DOMElement* findNextPmt(const DOMElement* rangeElt);

protected:
  /// A place to keep track of where we are within ACD data
  unsigned m_nDet;             // new
  unsigned m_nPmt;
  idents::AcdId  m_id;
  // don't think we need these
  //  unsigned m_nFace;
  //  unsigned m_nRow;
  //  unsigned m_nCol;
  //  unsigned m_nNA;

};

#endif
