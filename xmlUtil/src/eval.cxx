// $Header$
/*! \file Standalone program to transform source xml file into a 
    preprocessed version suitable for most clients (documentation-type
    clients will probably stick with the original file). 
  
    In particular this program will
      - add a <source> element to the output, indicating how and
        from what the output was generated
      - convert any <prim> values not already expressed in standard
        units (mm) to the standard
      - substitute values for all references constants, evaluating
        such constants as needed
      - delete the <derived> element and all descendents from the
        output
 */

#include "xml/XmlParser.h"
#include "xml/Dom.h"
#include "xmlUtil/Arith.h"
#include "xmlUtil/Substitute.h"
#include "xmlUtil/Constants.h"
#include "xmlUtil/Source.h"
#include "dom/DOM_Element.hpp"
#include "dom/DOM_NodeList.hpp"
#include "dom/DOM_DocumentType.hpp"

#include <string>
#include <iostream>
#include <fstream>

std::ostream *openOut(char * outfile);
void addSourceElt(DOM_Document doc);
void outProlog(const DOM_DocumentType& doctype, std::ostream& out);
char * stripDollar(char *toStrip);

const char chDoubleQ = 0x22;
const std::string dquote(&chDoubleQ);
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

  // prim processing
  xmlUtil::Constants *constants = new xmlUtil::Constants(doc);
  constants->normalizePrimary();

  // No longer need the evaluate step since substitution will
  // evaluate any constant which has not yet been evaluated.

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

  // Get rid of the whole <derived> part of <constants>.  No longer needed.
  DOM_Element derived = xml::Dom::findFirstChildByName(docElt, "constants");
  if (derived != DOM_Element() ) {
    derived = xml::Dom::findFirstChildByName(derived, "derived");

    if (derived != DOM_Element()) {
      xml::Dom::prune(derived);
      (derived.getParentNode()).removeChild(derived);
    }
  }

  // Add a <source> child to the outer gdd element
  xmlUtil::Source *source = 
    new xmlUtil::Source(doc, "xmlUtil/v1/src/eval.exe", "$Id$");
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




/*! Add a <source> element to the output gdd element to keep
 *  track of where it came from (CVS id of input) and how
 *  it was made (CVS id of this file)
 */
/*
void  addSourceElt(DOM_Document doc) {
  DOM_Element source = doc.createElement("source");
  DOM_NodeList gdds = doc.getElementsByTagName("gdd");
  // must be precisely one <gdd> element
  DOM_Node gddNode = gdds.item(0);
  DOM_Element& gdd = static_cast<DOM_Element &> (gddNode);

  DOMString idAtt = gdd.getAttribute("CVSid");
  const char * raw   = xml::Dom::transToChar(idAtt);
  unsigned   len = strlen(raw);
  char* myRaw = new char[len+1];
  strcpy(myRaw, raw);

  const  char *strippedRaw = stripDollar(myRaw);
  DOMString newAtt = DOMString(strippedRaw);

  source.setAttribute("inputId", newAtt);

  // And similarly for our own id.  Make a mutable copy.
  len = myId.size();
  char * creatorId = new char[len+1];
  strcpy(creatorId, myId.c_str());

  const char *creatorStripped = stripDollar(creatorId);
  source.setAttribute(DOMString("creatorId"), 
                      DOMString(creatorStripped));

  DOM_Text    sourceText = 
    doc.createTextNode("This file automatically created by xml/src/eval.cxx");

  source.appendChild(sourceText);

  delete [] creatorId;
  delete [] myRaw;
  // Finally need to insert this node as the first child of gdd
  gdd.insertBefore(source, gdd.getFirstChild());
}
*/

/*! Input is a null-terminated character string.  Output is
 * another such with leading and trailing '$', if any, stripped.
 * The purpose is to turn a CVS macro, such as "$Id$" into
 * something which will no longer be substituted for.
 */

/*
char * stripDollar(char* toStrip) {
  const char dollar = '$';

  char   nothing = 0;
  unsigned len = strlen(toStrip);

  if (dollar == toStrip[len - 1]) {
    len--;
    toStrip[len] = nothing;
  }
  if (dollar == toStrip[0]) {
    return (toStrip + 1);
  }
  return toStrip;
}
*/






