#include <xercesc/dom/DOM_Element.hpp>
#include <xercesc/dom/DOM_Document.hpp>
#include <xercesc/dom/DOMString.hpp>
#include "xml/Dom.h"
#include <fstream>


const char chDoubleQ = 0x22;
const std::string dquote(&chDoubleQ);

/*!
 *  Write out an xml declaration and copy the input DOCTYPE declaration
 *  to the output
 */
void outProlog(const DOM_DocumentType& doctype, std::ostream& out) {
  // write the xml declaration: <?xml version="1.0" ?>

  out << "<?xml version=" << dquote << "1.0" << dquote << "?>" << std::endl;
  if (doctype != DOM_DocumentType()) {
    char* transcoded = xml::Dom::transToChar(doctype.getName());
    if (transcoded != 0) {

      //    out << "<!DOCTYPE " << xml::Dom::transToChar(doctype.getName()) << " ";
      out << "<!DOCTYPE " << transcoded << " ";
    }
    else 
    {
      std::cout << "Failed to transcode doctype " << std::endl;
      return;
    }

    DOMString id = doctype.getPublicId();
    if (id != 0)   {
      //      out << " PUBLIC " << dquote << xml::Dom::transToChar(id) << dquote;
      transcoded = xml::Dom::transToChar(id);
      if (transcoded != 0) {
        out << " PUBLIC " << dquote << transcoded << dquote;
      }
      else  {
        std::cout << "Failed to transcode public id " << std::endl;
        return;
      }

      id = doctype.getSystemId();
      if (id != 0) {
        transcoded = xml::Dom::transToChar(id);
        if (transcoded != 0) {
          out << " " << dquote << transcoded << dquote;
        }
        else  {
          std::cout << "Failed to transcode system id " << std::endl;
          return;
        }
      }
    }
    else {
      id = doctype.getSystemId();
      if (id != 0)   {
        transcoded = xml::Dom::transToChar(id);
        if (transcoded != 0) {
          out << " SYSTEM " << dquote << transcoded << dquote;
        }
        else  {
          std::cout << "Failed to transcode system id " << std::endl;
          return;
        }
      }
    }
    id = doctype.getInternalSubset(); 
    if (id !=0) {
      transcoded = xml::Dom::transToChar(id);
      if (transcoded != 0) {
        //        out << "[" << xml::Dom::transToChar(id) << "]";
        out << "[" << transcoded << "]";
      }
      else {
        std::cout << "Failed to transcode internal subset" << std::endl;
        return;
      }
    }
    out << ">" << std::endl;
  }
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
