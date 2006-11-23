/** @file XmlFetchEvents.h

    @brief declaration of the FetchEvents class

$Header$

*/

#ifndef XmlFetchEvents_h
#define XmlFetchEvents_h


#include "IFetchEvents.h" 

#include <vector>
#include <xercesc/dom/DOMElement.hpp>

#include "xmlBase/XmlParser.h"
class TChain;
class TTree;
class TFile;


/** @class XmlFetchEvents
    @brief manage the retrieval of events using some specified parameter(s) from an XML file
    @author Heather Kelly heather625@gmail.com


*/
class XmlFetchEvents : public IFetchEvents
{
public:

    XmlFetchEvents(const std::string& xmlFile, const std::string& param);

    virtual ~XmlFetchEvents();

    double getAttributeValue(const std::string& elemName, double binVal);

    int getFiles(double binVal, TChain* chain);

    /// Purpose and Method:  Returns a TTree constructed from the "fileList" associated with the bin
    /// found using binVal.
    /// Returns a pointer to the TTree  if completely successful
    /// returns 0 if still in the same bin as a previous call
    TTree* getTree(double binVal);

    virtual double minVal()const{return m_minval;}///< return minimum value allowed
    virtual double maxVal()const{return m_maxval;}///< return maximum value allowed

    static double m_badVal;

private:

    /// actually handles the XML reading
    xmlBase::XmlParser m_xmlParser;

    /// stores all children of the top element in the XML file
    std::vector<xmlBase::DOMElement* > m_children;  

    /// children that contain the parameter we're interested in
    std::vector<xmlBase::DOMElement*> m_paramChildren; 

    /// children that contain the bins for the parameter we're interested in
    std::vector<xmlBase::DOMElement*> m_binChildren;

    /// Store the most recently accessed DOMElement in the m_binChildren vector;
    int m_lastBinIndex;
    double m_lastBinMin, m_lastBinMax;

    double m_minval, m_maxval;

    TFile * m_file; ///< current TFile
};


#endif
