#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMDocumentType.hpp>
#include <xercesc/util/XMLString.hpp>
#include <iostream>
#include <fstream>


const char chDoubleQ = 0x22;
const std::string dquote(&chDoubleQ);

/*!
 *  Write out an xml declaration and copy the input DOCTYPE declaration
 *  to the output
 */
void outProlog(const XERCES_CPP_NAMESPACE_QUALIFIER DOMDocumentType* doctype, 
               std::ostream& out) {
  // write the xml declaration: <?xml version="1.0" ?>
  using XERCES_CPP_NAMESPACE_QUALIFIER XMLString;

  out << "<?xml version=" << dquote << "1.0" << dquote << "?>" << std::endl;
  if (doctype != 0) {
    char* transcoded = XMLString::transcode(doctype->getName());
    if (transcoded != 0) {

      out << "<!DOCTYPE " << transcoded << " ";
      XMLString::release(&transcoded);
    }
    else 
    {
      std::cout << "Failed to transcode doctype " << std::endl;
      return;
    }

    const XMLCh* id = doctype->getPublicId();
    if (id != 0)   {

      transcoded = XMLString::transcode(id);
      if (transcoded != 0) {
        out << " PUBLIC " << dquote << transcoded << dquote;
        XMLString::release(&transcoded);
      }
      else  {
        std::cout << "Failed to transcode public id " << std::endl;
        return;
      }

      id = doctype->getSystemId();
      if (id != 0) {
        transcoded = XMLString::transcode(id);
        if (transcoded != 0) {
          out << " " << dquote << transcoded << dquote;
          XMLString::release(&transcoded);
        }
        else  {
          std::cout << "Failed to transcode system id " << std::endl;
          return;
        }
      }
    }
    else {
      id = doctype->getSystemId();
      if (id != 0)   {
        transcoded = XMLString::transcode(id);
        if (transcoded != 0) {
          out << " SYSTEM " << dquote << transcoded << dquote;
          XMLString::release(&transcoded);
        }
        else  {
          std::cout << "Failed to transcode system id " << std::endl;
          return;
        }
      }
    }
    id = doctype->getInternalSubset(); 
    if (id !=0) {
      transcoded = XMLString::transcode(id);
      if (transcoded != 0) {
        out << "[" << transcoded << "]";
        XMLString::release(&transcoded);
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
