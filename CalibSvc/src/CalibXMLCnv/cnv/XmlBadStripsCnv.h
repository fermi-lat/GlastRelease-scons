// $Header$
#ifndef CalibSvc_XmlBadStripsCnv_h
#define CalibSvc_XmlBadStripsCnv_h

/** @class XmlBadStripsCnv_h

  Converter from xml for tracker bad strips: collections declared hot or dead

  @author J. Bogart
*/
#include "XmlBaseCnv.h"
#include <dom/DOM_Element.hpp>
#include "CalibData/Tkr/BadStrips.h"


template <class TYPE> class CnvFactory;

class XmlBadStripsCnv : public XmlBaseCnv {

  /// Friend needed for instantiation
  friend class CnvFactory<XmlBadStripsCnv>;
public:
  const CLID& objType() const;
  static const CLID& classID();
protected:

  XmlBadStripsCnv(ISvcLocator* svcs);

  virtual ~XmlBadStripsCnv() {}

  /** Does the real work of conversion.
      @param  element     Document element of XML document (input)
      @param  refpObject  Pointer to newly-created object goes here (output)
  */
  virtual StatusCode i_createObj(const DOM_Element& element,
                                 DataObject*& refpObject);


  /// Internal helper; handles xml for a tower
  StatusCode processTower(const DOM_Element& towerElt, 
                          CalibData::BadStrips* pBad);

  /** Internal helper; handles xml for a uniplane
      @param uniElt  DOM element representing bad strips within a (uni)plane
      @param row     row of tower containing the uniplane (input)
      @param col     column of tower containing the uniplane (input)
      @param pBad    Pointer to BadStrips object to be updated
   */
  StatusCode processUni(const DOM_Element& uniElt,
                        unsigned row, unsigned col,
                        CalibData::BadStrips* pBad);

  StatusCode addStrips(const DOM_Element& uniElt, CalibData::StripCol* strips);


};
#endif
