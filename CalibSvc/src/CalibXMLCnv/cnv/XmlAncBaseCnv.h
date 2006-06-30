// $Header$
#ifndef CalibData_XmlAncBaseCnv_h
#define CalibData_XmlAncBaseCnv_h

/** @class XmlAncBaseCnv 

  Base class for ANC calibration converters from XML files to TCDS.
  All such converters need to do certain things, which are
  handled here.  Methods common to *all* calibrations are in the
  base class XmlBaseCnv

  @author J. Bogart
*/

#include "XmlBaseCnv.h"
class XmlAncBaseCnv : public XmlBaseCnv {

public:

  XmlAncBaseCnv(ISvcLocator* svc, const CLID& clid);

  virtual ~XmlAncBaseCnv() {};


  /** Convenience routine used by most ANC calibration types, which
   have a <dimension> element describing how the remainder of the
   data is laid out.
   @a nMod  number of modules
   @a nLay  number of layers per module
   @a nChan number of channels per layer
       Note "layer" really only has meaning for tagger.  qdc modules
       don't have layers, but treat as if nLay is
  */
  StatusCode readAncDimension(const DOMElement* docElt, unsigned& nMod,
                              unsigned& nLay, unsigned& nChan);

  /// Another one to find first ped, gain, etc. (depending on calib. type)
  DOMElement* findFirstChan(const DOMElement* docElt);

  /// Still another one to navigate XML file and find next set of pmt data
  DOMElement* findNextChan(const DOMElement* rangeElt);

protected:
  /// A place to keep track of where we are within ANC data
  unsigned m_iMod;
  unsigned m_iLay;
  unsigned m_iChan;
  /// Tagger files have <layer> element; qdc files don't
  bool  m_hasLayers;
  bool  m_isTagger;
  bool  m_isQdc;
  std::string m_eltName;
};

#endif
