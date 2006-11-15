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
//#include "xmlBase/Dom.h"
//#include <xercesc/dom/DOMElement.hpp>
//#include <xercesc/dom/DOMNodeList.hpp>


/** @class XmlFetchEvents
    @brief manage the retrieval of events using some specified parameter(s)
    @author Heather Kelly heather625@gmail.com


*/
class XmlFetchEvents : public IFetchEvents
{
public:

    /** @brief ctor
        @param param parameter used for retrieving data from xmlFile
        @param xmlFile name of the file containing the event store
    */
    //XmlFetchEvents(const std::string& xmlFile);

    XmlFetchEvents(const std::string& xmlFile, const std::string& param);

    ~XmlFetchEvents();

    double getAttributeValue(const std::string& elemName, double binVal);

    int getFiles(double binVal, TChain* chain);


private:

  //  std::string m_param; /// name of the variable that we use to key the rates

    xmlBase::XmlParser m_xmlParser;

    std::vector<xmlBase::DOMElement* > m_children;  // stores all children of top element in XML file

    std::vector<xmlBase::DOMElement*> m_paramChildren; // children that contain the parameter we're interested in

    std::vector<xmlBase::DOMElement*> m_binChildren;

    /// Store the most recently accessed DOMElement in the m_binChildren vector;
    int m_lastBinIndex;
    double m_lastBinMin, m_lastBinMax;


};


#endif
