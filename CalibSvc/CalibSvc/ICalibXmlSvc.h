// $Header$

#ifndef IXMLSvc_h 
#define IXMLSvc_h

/** @class ICalibXmlSvc
   Abstract interface to be satisfied by any XML conversion implementation.
 
   Will parse a file, making its DOM rep. available

   Maybe also provide layer over "raw" access to the DOM to insulate
   converters from changes moving from one xml parser version (or even
   parser) to another.
   ...or do we just depend on things in calibUtil?

*/
#include "GaudiKernel/IInterface.h"
#include <dom/DOM_Document.hpp>

static const InterfaceID IID_ICalibXmlSvc("ICalibXmlSvc", 1, 0);

class ICalibXmlSvc : virtual public IInterface 
{
public:
  // Re-implemented from IInterface
  static const InterfaceID& interfaceID() { return IID_ICalibXmlSvc; }
  
  /**
   * This method parses an xml file and produces the corresponding DOM
   * document.
   * @param fileName the name of the file to parse
   * @return the document issued from the parsing
   */
  virtual DOM_Document parse(const char* filename) = 0;

  // Do we also want a "reset" or "clearDocument" ?  Can in any case
  // do this internally when a new document is to be parsed so might not
  // be necessary to have explicit public method.
};


#endif
