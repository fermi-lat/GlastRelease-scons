// $Header$
/*! \file Standalone program to transform source xml file into a 
    preprocessed version suitable for documentation.  Typically
    will be transformed further, e.g. to html by an xslt transform.
    xslt is not suitable to do the whole job because it can't handle
    recursive evaluation of the constants.
  
    In particular this program will
      - add a <source> element to the output, indicating how and
        from what the output was generated
      - evaluate all <const>
      - delete all content except <notes> elements from each <const>
      - delete all <section> elements (for now.  Later may think of
        some useful information for doc. purposes to be extracted
        from them).
      - output the resulting xml document.
 */

#include "xml/XmlParser.h"
#include "xml/Dom.h"
#include "xmlUtil/Source.h"
#include "xmlUtil/Arith.h"
#include "xmlUtil/Constants.h"
#include "dom/DOM_Element.hpp"
#include "dom/DOM_NodeList.hpp"
#include "dom/DOM_DocumentType.hpp"

#include <string>
#include <iostream>
#include <fstream>

std::ostream *openOut(char * outfile);

void outProlog(const DOM_DocumentType& doctype, std::ostream& out);

char * stripDollar(char *toStrip);

const char chDoubleQuote = 0x22;
const std::string dquote(&chDoubleQuote);
const std::string myId("$Id$");

/*!
    Main program for the forDoc application.
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

  // For each const child of a derCategory
  //   evaluate it, saving value
  //   delete all children except <notes>
  xmlUtil::Constants *constants = new xmlUtil::Constants(doc);
  constants->evalConstants();
  constants->pruneConstants(true);  

  // Delete all sections
  DOM_NodeList sections = docElt.getElementsByTagName(DOMString("section"));
  DOM_Node secNode = sections.item(0);

  while (secNode != DOM_Node() ) {
    DOM_Node toCome = secNode.getNextSibling();
    DOM_Element& secElt = static_cast<DOM_Element &> (secNode);
    xml::Dom::prune(secElt);
    (secElt.getParentNode()).removeChild(secElt);
    secNode = toCome;
  }


  // Add a <source> child to the outer gdd element
  xmlUtil::Source * source = 
    new xmlUtil::Source(doc, "xmlUtil/v1/src/forDoc.exe", myId.c_str());
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
