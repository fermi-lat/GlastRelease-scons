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

class XmlAcdBaseCnv : public XmlBaseCnv {

public:

  XmlAcdBaseCnv(ISvcLocator* svc, const CLID& clid);

  virtual ~XmlAcdBaseCnv() {};


  /// Convenience routine used by most ACD calibration types, which
  /// have a <dimension> element describing how the remainder of the
  /// data is laid out.
  StatusCode readAcdDimension(const DOMElement* docElt,
                              unsigned& nFace, unsigned& nRow,
                              unsigned& nCol, unsigned& nPmt,
                              unsigned& nRange);


  /// Another one to find first range element
  DOMElement* findFirstRange(const DOMElement* docElt);

  /// Still another one to navigate XML file and find next set of range data
  DOMElement* findNextRange(const DOMElement* rangeElt);

protected:
  /// A place to keep track of where we are within ACD data
  unsigned m_nFace;
  unsigned m_nRow;
  unsigned m_nCol;
  unsigned m_nPmt;
  unsigned m_nRange;

};

#endif
