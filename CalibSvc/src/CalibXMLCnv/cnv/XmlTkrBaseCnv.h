// $Header$
#ifndef CalibData_XmlTkrBaseCnv_h
#define CalibData_XmlTkrBaseCnv_h

/** @class XmlTkrBaseCnv 

  Base class for *most* TKR calibration converters from XML files to TCDS.
  (bad strips are an exception; converters do not inherit from this class).
  All such converters need to do certain things, which are
  handled here.  Methods common to *all* calibrations are in the
  base class XmlBaseCnv

  @author J. Bogart
*/

#include "XmlBaseCnv.h"

class XmlTkrBaseCnv : public XmlBaseCnv {

public:

  XmlTkrBaseCnv(ISvcLocator* svc, const CLID& clid);

  virtual ~XmlTkrBaseCnv() {};


  /// Convenience routine used by most CAL calibration types, which
  /// have a <dimension> element describing how the remainder of the
  /// data is laid out.
  StatusCode readDimension(const DOM_Element& docElt, 
                           unsigned& nRow, unsigned& nCol, 
                           unsigned& nTray,
                           unsigned& nChip);

  unsigned countChannels(const DOM_Element& docElt, std::string eltName);
  // Another one to find first range element
  //  DOM_Element findFirstRange(const DOM_Element& docElt);

  // Still another one to navigate XML file and find next set of range data
  //DOM_Element findNextRange(const DOM_Element& rangeElt);

protected:
  // Number of channels having data in persistent rep.
  unsigned m_nChannel;

  // A place to keep track of where we are if we're handling CAL data
  /*
  unsigned m_nRow;
  unsigned m_nCol;
  unsigned m_nLayer;
  unsigned m_nXtal;
  unsigned m_nFace;
  unsigned m_nRange;
  */

};

#endif
