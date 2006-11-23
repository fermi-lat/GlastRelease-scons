/**  @file XmlFetchEvents.cxx
@brief implementation of class XmlFetchEvents

$Header$  
*/

#include "XmlFetchEvents.h"
#include "xmlBase/Dom.h"
#include "facilities/Util.h"
#include <xercesc/dom/DOMNodeList.hpp>

#include "TChain.h"
#include "TFile.h"

#include <stdexcept>
#include <sstream>
#include <cmath>
#include <cassert>


using XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument;
using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;
using xmlBase::Dom;

double XmlFetchEvents::m_badVal = -999999.0;

XmlFetchEvents::XmlFetchEvents(const std::string& xmlFile, const std::string& param)
: IFetchEvents(xmlFile,param) 
, m_file(0)
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

    // retrieve the top-level element in our XML 
    DOMElement* elem = doc->getDocumentElement();
    //xmlBase::Dom::prettyPrintElement(elem, std::cout, "");

    // Assuming "source" is the most basic element in our XML file
    xmlBase::Dom::getChildrenByTagName(elem, "source", m_children);
    m_paramChildren.clear();
    m_paramChildren.reserve(m_children.size()); // reserve enough space for all possible children

    // Save all children elements that pertain to our search parameter
    std::vector<DOMElement*>::const_iterator domElemIt; 
    for (domElemIt = m_children.begin(); domElemIt != m_children.end(); domElemIt++) { 
        std::string typeAttr = xmlBase::Dom::getAttribute(*domElemIt, "type");
        if (typeAttr == m_param) {
            // found one save it
            m_paramChildren.push_back(*domElemIt);
            double minval(facilities::Util::stringToDouble(xmlBase::Dom::getAttribute(*domElemIt, "rangeMin"))),
                maxval(facilities::Util::stringToDouble(xmlBase::Dom::getAttribute(*domElemIt, "rangeMax")));
            if(minval< m_minval) m_minval = minval;
            if(maxval>m_maxval) m_maxval = maxval;
        }
    }

    // Save all the "bin" elements associated with our parameter
    m_binChildren.clear();
    for(domElemIt=m_paramChildren.begin(); domElemIt != m_paramChildren.end(); domElemIt++) {
        // tell method not to clear our binChildren vector by setting last param to false
        xmlBase::Dom::getChildrenByTagName(*domElemIt, "bin", m_binChildren, false);
    }

    if( m_binChildren.empty() ){
        throw std::invalid_argument("XmlFetchEvents: did not find entries for "+param);
    }
    m_lastBinIndex = -1;
    m_lastBinMin = m_badVal;
    m_lastBinMax = m_badVal;
}


XmlFetchEvents::~XmlFetchEvents()
{
    m_children.clear();
    m_paramChildren.clear();
    m_binChildren.clear();
    delete m_file;
}

double XmlFetchEvents::getAttributeValue(const std::string& elemName, double binVal) {
    /// Purpose and Method:  Extracts any attribute associated with the name elemName from
    /// our vector of bins.  Uses the binVal to determine which bin to use.
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
    
    return m_badVal;
}


int XmlFetchEvents::getFiles(double binVal, TChain* chain) {
    /// Purpose and Method:  Returns a TChain constructed from the "fileList" associated with the bin
    /// found using binVal.
    /// Returns 0 if completely successful
    /// Returns 1 if any of files failed to be added to the TChain
    /// Returns -1 if no files were found and the TChain remains empty

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
            int status = (dynamic_cast<TChain*>(chain))->AddFile(fileNameStr.c_str(),0,treeNameStr.c_str());
            if (status == 0) statFlag |= 1;
        }
    } else
        statFlag = -1;

    return(statFlag);
}

TTree* XmlFetchEvents::getTree(double binVal) {
    /// Purpose and Method:  Returns a TTree constructed from the "fileList" associated with the bin
    /// found using binVal.
    /// Returns a pointer to the TTree  if  successful
    /// Return 0 if the TTree was previously found
    /// throws an exception if failure

    // Check to see if we are still in the same bin that was previously set up

    if (m_lastBinIndex >= 0 && (binVal >= m_lastBinMin) && (binVal < m_lastBinMax) ) {
        std::cout << "no new tree needed" << std::endl;
        return 0;  // no new ttree to find.
    }
    
    
    std::vector<DOMElement*> fileList;
    TTree* tree(0);
  
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
            break;
        }
        ++curIndex;
    }
    if( fileList.empty() ){
        std::stringstream err;
        err << "XMLFetchEvents::getTree -- no files found for "<< m_param << "=" << binVal;

        throw std::runtime_error(err.str());
    }

    if (fileList.size() ==1) {
        DOMElement* file= fileList.front();
        std::string fileNameStr = xmlBase::Dom::getAttribute(file, "filePath");
        facilities::Util::expandEnvVar(&fileNameStr);
        std::string treeNameStr = xmlBase::Dom::getAttribute(file, "treeName");
        delete m_file;  // get rid of last guy
        m_file = new TFile(fileNameStr.c_str());
        tree = (TTree*)m_file->Get(treeNameStr.c_str());
        if( tree==0 ) throw std::runtime_error("XMLFetchEvent: did not find the TTree "+treeNameStr);
    }else{
        throw std::runtime_error("XMLFetchEvents::getTree expected a single file");
    }
    return tree;
}
