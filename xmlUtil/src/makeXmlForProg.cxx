// $Header$
/*! \file Standalone program to transform source xml file into a 
    preprocessed version suitable for most clients programs (such
    as Simulation and Reconstruction).     Clients needing to
    extract documentation from the xml file will go a different
    route (see makeXmlForDoc.cxx.)

  
    In particular this program will
      - add a <source> element to the output, indicating how and
        from what the output was generated
      - convert any <prim> values not already expressed in standard
        units (mm) to the standard
      - substitute values for all references constants, evaluating
        such constants as needed
      - delete any <derCategory> without the "save" attribute set
        to true and all its children.  For any <derCategory> which
        is to be saved, delete all element content from each
        <const> (except for notes??).
 */

#include "xml/XmlParser.h"
#include "xml/Dom.h"
#include "xmlUtil/Source.h"
#include "xmlUtil/Constants.h"
#include "xmlUtil/Arith.h"
#include "xmlUtil/Substitute.h"
#include "dom/DOM_Element.hpp"
#include "dom/DOM_NodeList.hpp"
#include "dom/DOM_NamedNodeMap.hpp"
#include "dom/DOM_DocumentType.hpp"

#include <string>
#include <iostream>
#include <fstream>

std::ostream *openOut(char * outfile);
// void addSourceElt(DOM_Document doc);
void outProlog(const DOM_DocumentType& doctype, std::ostream& out);
// char * stripDollar(char *toStrip);

const char chDoubleQuote = 0x22;
const std::string dquote(&chDoubleQuote);
const std::string myId("$Id$");

/*!
    Main program for the eval application.
    \param arg1 is the input xml file
    \param arg2 is specification for the output file, or "-" to output
    to std::cout
*/
int main(int argc, char* argv[]) {
  std::ostream *out;
  if (argc < 3) {  // instructions
    std::cout << "Required first argument is xml file to be parsed" 
              << std::endl;
    std::cout << "Required second argument is output file (- for stdout)"
              << std::endl;
    exit(0);
  }

  xml::XmlParser* parser = new xml::XmlParser();
  DOM_Document doc = parser->parse(argv[1]);

  if (doc == 0) {
    std::cout << "Document failed to parse correctly" << std::endl;
    exit(0);
  }

  // Else successful.   Open output
  out = openOut(argv[2]);

  DOM_Element docElt = doc.getDocumentElement();
  DOM_DocumentType  doctype = doc.getDoctype(); //  check for gdd?? 
  xmlUtil::Constants   *constants = new xmlUtil::Constants(doc);

  // prim processing
  constants->normalizePrimary();

  // Evaluate all constants.  Most would be evaluated even without
  // this step because of the use of Substitute, but if there happens
  // to be a constant not actually used in the <section> part of
  // the file, Substitute would not cause that constant to be evaluated.
  constants->evalConstants();

  // Search for all dimension and position references, substitute
  // value of referenced element.
  DOM_NodeList sections = docElt.getElementsByTagName(DOMString("section"));
  int nSec = sections.getLength();
  int iSec;

  xmlUtil::Substitute* sub = new xmlUtil::Substitute(doc);

  for (iSec = 0; iSec < nSec; iSec++) {
    int nSub;
    DOM_Node  secNode = sections.item(iSec);
    DOM_Element& secElt = static_cast<DOM_Element &> (secNode);
    nSub = sub->execute(secElt);
    std::cout << "#elements substituted for in this section:  " << nSub
              << std::endl;
  }
  
  // Get rid of any <derCategory> without the save attribute set to "true".
  // For the rest, prune element content from all <const> children;
  // it's no longer needed.
  DOM_Element derived = xml::Dom::findFirstChildByName(docElt, "constants");
  if (derived != DOM_Element() ) {
    derived = xml::Dom::findFirstChildByName(derived, "derived");
    
    if (derived != DOM_Element()) {
      DOM_NodeList cats = derived.getElementsByTagName("derCategory");
      int nCat = cats.getLength();
      int iCat;
      for (iCat=0; iCat < nCat; iCat++) {
        DOM_Node node = cats.item(iCat);
        DOM_Node attNode = (node.getAttributes()).getNamedItem("save");
        if ((attNode.getNodeValue()).equals(DOMString("false")) ) {  // remove
          DOM_Element& toPrune = static_cast<DOM_Element &>(node);
          xml::Dom::prune(toPrune);
          (toPrune.getParentNode()).removeChild(toPrune);
          // compensate for changing node list
          iCat--;  nCat--;
        }
      }

      // Now just get rid of all content for <const> children 
      // in remaining derived categories
      constants->pruneConstants();
    }
  }

  // Add a <source> child to the outer gdd element
  xmlUtil::Source * source = 
    new xmlUtil::Source(doc, "xmlUtil/v1/src/forProg.exe", myId.c_str());
  source->add();
  
  // Output the xml declaration and all the text in the DOCTYPE (see DOMPrint)
  outProlog(doctype, *out);

  // Finally output the elements
  // May want option to exclude comments here
  xml::Dom::prettyPrintElement(docElt, *out, "");

  delete parser;
  return(0);
}

/*! 
 *  Open specified output file (may be standard output if "-"
 *  is given as input argument)
 */
std::ostream *openOut(char * outfile)
{
  char  *hyphen = "-";
  std::ostream* out;
  
  if (*outfile == *hyphen) {
    out = &std::cout;
  }
  else {   // try to open file as ostream
    out = new std::ofstream(outfile);
  }
  return out;
}

/*!
 *  Write out an xml declaration and copy the input DOCTYPE declaration
 *  to the output
 */
void outProlog(const DOM_DocumentType& doctype, std::ostream& out) {
  // write the xml declaration: <?xml version="1.0" ?>

  out << "<?xml version=" << dquote << "1.0" << dquote << "?>" << std::endl;
  if (doctype != DOM_DocumentType()) {

    out << "<!DOCTYPE " << xml::Dom::transToChar(doctype.getName()) << " ";

    DOMString id = doctype.getPublicId();
    if (id != 0)   {
      out << " PUBLIC " << dquote << xml::Dom::transToChar(id) << dquote;
      id = doctype.getSystemId();
      if (id != 0) {
        out << " " << dquote << xml::Dom::transToChar(id) << dquote;
      }
    }
    else {
      id = doctype.getSystemId();
      if (id != 0)   {
        out << " SYSTEM " << dquote << xml::Dom::transToChar(id) << dquote;
      }
    }
    id = doctype.getInternalSubset(); 
    if (id !=0) {
      out << "[" << xml::Dom::transToChar(id) << "]";
    }
    out << ">" << std::endl;
  }
}
