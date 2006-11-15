/**  @file FetchEvents.cxx
    @brief implementation of class FetchEvents
    
  $Header$  
*/

#include "XmlFetchEvents.h"
#include "xmlBase/Dom.h"
#include "facilities/Util.h"
#include <xercesc/dom/DOMNodeList.hpp>

    using XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument;
    using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;
    using xmlBase::Dom;

    /*
XmlFetchEvents::XmlFetchEvents(const std::string& xmlFile)
: IFetchEvents(xmlFile)
{
   XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument* doc;
   m_xmlParser.doSchema(true);

    try {
        doc = m_xmlParser.parse(m_dataStore.c_str());
    } catch (std::exception ex) {
        std::cerr << ex.what();
        std::cerr.flush();
        throw ex;
    }
    if (doc == NULL) {
      throw std::runtime_error("Archive does not have proper myHost.xml file");
    }

   DOMElement* elem = doc->getDocumentElement();
   xmlBase::Dom::prettyPrintElement(elem, std::cout, "");
   
   xmlBase::Dom::getChildrenByTagName(elem, "source", m_children);

   m_lastBinIndex = -1;
   m_lastBinMin = -99999.;
   m_lastBinMax = -99999.;
}*/

XmlFetchEvents::XmlFetchEvents(const std::string& xmlFile, const std::string& param) : IFetchEvents(xmlFile,param) {
    //m_param = param;
   XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument* doc;
   m_xmlParser.doSchema(true);

    try {
        doc = m_xmlParser.parse(m_dataStore.c_str());
    } catch (std::exception ex) {
        std::cerr << ex.what();
        std::cerr.flush();
        throw ex;
    }
    if (doc == NULL) {
      throw std::runtime_error("Archive does not have proper myHost.xml file");
    }

   DOMElement* elem = doc->getDocumentElement();
   xmlBase::Dom::prettyPrintElement(elem, std::cout, "");
   
   xmlBase::Dom::getChildrenByTagName(elem, "source", m_children);
   m_paramChildren.clear();
   m_paramChildren.reserve(m_children.size()); // reserve enough space for all possible children

   // Save all children elements that pertain to our search parameter
    std::vector<DOMElement*>::const_iterator domElemIt; 
    for (domElemIt = m_children.begin(); domElemIt != m_children.end(); domElemIt++) { 
       std::string typeAttr = xmlBase::Dom::getAttribute(*domElemIt, "type");
       if (typeAttr == m_param) {
           m_paramChildren.push_back(*domElemIt);
       }
   }

   m_binChildren.clear();
   for(domElemIt=m_paramChildren.begin(); domElemIt != m_paramChildren.end(); domElemIt++) {
       // tell method not to clear our binChildren vector by setting last param to false
       xmlBase::Dom::getChildrenByTagName(*domElemIt, "bin", m_binChildren, false);
   }

   m_lastBinIndex = -1;
   m_lastBinMin = -99999.;
   m_lastBinMax = -99999.;
}


XmlFetchEvents::~XmlFetchEvents()
{
    m_children.clear();
    m_paramChildren.clear();
    m_binChildren.clear();
}

double XmlFetchEvents::getAttributeValue(const std::string& elemName, double binVal) {
    // Check to see if we're accessing the same bin we have previously
    if (m_lastBinIndex >= 0) {
       if ((binVal >= m_lastBinMin) && (binVal <= m_lastBinMax)) {
           std::string retValStr = xmlBase::Dom::getAttribute(m_binChildren[m_lastBinIndex], elemName);
           double retVal = facilities::Util::stringToDouble(retValStr);
           return retVal;
       }
    }
    // Otherwise search the whole vector
    std::vector<DOMElement*>::const_iterator domElemIt; 
    int curIndex = 0;

    for (domElemIt = m_binChildren.begin(); domElemIt != m_binChildren.end(); domElemIt++) { 
       std::string minBinStr = xmlBase::Dom::getAttribute(*domElemIt, "min");
       std::string maxBinStr = xmlBase::Dom::getAttribute(*domElemIt, "max");
       double minBin = facilities::Util::stringToDouble(minBinStr);
       double maxBin = facilities::Util::stringToDouble(maxBinStr);
       if ((binVal >= minBin) && (binVal <= maxBin)) {
           m_lastBinIndex = curIndex;
           m_lastBinMin = minBin;
           m_lastBinMax = maxBin;
           std::string retValStr = xmlBase::Dom::getAttribute(*domElemIt, elemName);
           double retVal = facilities::Util::stringToDouble(retValStr);
           return retVal;
       }
       ++curIndex;
   }
   return -9999999.0;
}


int XmlFetchEvents::getFiles(double binVal, TChain* chain) {
    int statFlag = 0;
    std::vector<DOMElement*> fileList;
    fileList.clear();
    // Check to see if we're accessing the same bin we have previously
    if (m_lastBinIndex >= 0) {
       if ((binVal >= m_lastBinMin) && (binVal <= m_lastBinMax)) {
           DOMElement* fileListElem = xmlBase::Dom::findFirstChildByName(m_binChildren[m_lastBinIndex], "fileList");
           xmlBase::Dom::getChildrenByTagName(fileListElem, "file", fileList);
       }
    } else {
    // Otherwise search the whole vector
    std::vector<DOMElement*>::const_iterator domElemIt; 
    int curIndex = 0;

    for (domElemIt = m_binChildren.begin(); domElemIt != m_binChildren.end(); domElemIt++) { 
       std::string minBinStr = xmlBase::Dom::getAttribute(*domElemIt, "min");
       std::string maxBinStr = xmlBase::Dom::getAttribute(*domElemIt, "max");
       double minBin = facilities::Util::stringToDouble(minBinStr);
       double maxBin = facilities::Util::stringToDouble(maxBinStr);
       if ((binVal >= minBin) && (binVal <= maxBin)) {
           m_lastBinIndex = curIndex;
           m_lastBinMin = minBin;
           m_lastBinMax = maxBin;
           DOMElement* fileListElem = xmlBase::Dom::findFirstChildByName(*domElemIt, "fileList");
           xmlBase::Dom::getChildrenByTagName(fileListElem, "file", fileList);
       }
       ++curIndex;
   }
    }
  
    if (fileList.size() > 0) {
        std::vector<DOMElement*>::const_iterator fileIt;
        for (fileIt=fileList.begin(); fileIt != fileList.end(); fileIt++) {
            std::string fileNameStr = xmlBase::Dom::getAttribute(*fileIt, "filePath");
            facilities::Util::expandEnvVar(&fileNameStr);
            std::string treeNameStr = xmlBase::Dom::getAttribute(*fileIt, "treeName");
            // Returns 1 if success and 0 if there is a problem adding the file
            statFlag = chain->AddFile(fileNameStr.c_str(),0,treeNameStr.c_str());
        }
    }

   return(statFlag);
}
