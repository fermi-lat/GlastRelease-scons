/**@file xmlSplitEngineFactory.cxx

@brief implementation of class xmlSplitEngineFactory

$Header$
*/

#include "xmlSplitEngineFactory.h"
#include "../ImActivityNodes/SplitEngineNode.h"
#include "xmlBase/XmlParser.h"
#include "facilities/Util.h"
#include <xercesc/dom/DOMElement.hpp>
#include "xmlBase/Dom.h"

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include <set> 
#include <cmath>  // for M_PI, among others

XERCES_CPP_NAMESPACE_USE


namespace {

    // output errors to console. (note: m_ prefex for data members )
    // only for statistics to show in debug output
    int  nesting_level, max_depth, node_count;  
    double min_prob, max_prob;
} // anonomous namespace

xmlSplitEngineFactory::xmlSplitEngineFactory(std::ostream& log, int iVerbosity )
                        : xmlFactoryBase(log,iVerbosity), 
                          m_log(log), m_outputLevel(iVerbosity)
{
}

IImActivityNode* xmlSplitEngineFactory::operator()(const DOMElement* xmlActivityNode)
{
    // Retrieve name and node id
    DOMElement* displayInfo = xmlBase::Dom::findFirstChildByName(xmlActivityNode, "DisplayInfo");
    std::string sType       = "SplitEngineNode";
    std::string sName       = xmlBase::Dom::getAttribute(displayInfo, "labelText");
    std::string sId         = xmlBase::Dom::getAttribute(xmlActivityNode, "id");

    // Create the node
    SplitEngineNode* node = new SplitEngineNode(sType, sName, sId);

    DOMElement* xmlProperty = getXTProperty(xmlActivityNode, "testExpression");

    std::string sExpression;

    try
    {
        DOMElement* xmlComplex = xmlBase::Dom::findFirstChildByName(xmlProperty, "Complex");
        sExpression = xmlBase::Dom::getTextContent(xmlComplex);
    }
    //catch(xmlBase::WrongNodeType& e)
    catch(...)
    {
        sExpression = xmlBase::Dom::getAttribute(xmlProperty, "value");
    }

    // Trim the blank spaces
    sExpression = trimBlanks(sExpression);
        
    // Store
    node->setExpression(sExpression);

    // Parse
    StringList parsedExpression;
    parseExpression(parsedExpression, sExpression);

    node->setParsedExpression(parsedExpression);

    return node;
}

xmlSplitEngineFactory::~xmlSplitEngineFactory()
{
    ///@TODO: delete trees
}
