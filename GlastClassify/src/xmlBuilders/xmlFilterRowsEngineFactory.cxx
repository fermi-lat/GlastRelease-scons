/**@file xmlFilterRowsEngineFactory.cxx

@brief implementation of class xmlFilterRowsEngineFactory

$Header$
*/

#include "xmlFilterRowsEngineFactory.h"
#include "../ImActivityNodes/FilterRowsEngineNode.h"
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

xmlFilterRowsEngineFactory::xmlFilterRowsEngineFactory(XTExprsnParser& parser) : xmlFactoryBase(parser)
{
}

IImActivityNode* xmlFilterRowsEngineFactory::operator()(const DOMElement* xmlActivityNode)
{
    // Retrieve name and node id
    DOMElement* displayInfo = xmlBase::Dom::findFirstChildByName(xmlActivityNode, "DisplayInfo");
    std::string sType       = "FilterRowsEngineNode";
    std::string sName       = xmlBase::Dom::getAttribute(displayInfo, "labelText");
    std::string sId         = xmlBase::Dom::getAttribute(xmlActivityNode, "id");

    // Create the node
    FilterRowsEngineNode* node = new FilterRowsEngineNode(sType, sName, sId);

    // Decode the expression
    DOMElement* xmlProperty = getXTProperty(xmlActivityNode, "testExpression");
    std::string sExpression;

    try
    {
        DOMElement* xmlComplex  = xmlBase::Dom::findFirstChildByName(xmlProperty, "Complex");
        sExpression = xmlBase::Dom::getTextContent(xmlComplex);
    }
    //catch(xmlBase::WrongNodeType& e)
    catch(...)
    {
        sExpression = xmlBase::Dom::getAttribute(xmlProperty, "value");
    }

    // Store
    node->setExpression(sExpression);

    // Parse
    IXTExprsnNode* xprsnNode = XprsnParser().parseExpression(sExpression);

    // Save the executable parsed expression
    node->setXTExprsnNode(xprsnNode);

    // Ok, now determine if we are "including" or "excluding" 
    xmlProperty = getXTProperty(xmlActivityNode, "includeRows");

    sExpression = xmlBase::Dom::getAttribute(xmlProperty, "value");

    bool keepRow = sExpression == "true";
    node->setIncludeRows(keepRow);

    return node;
}

xmlFilterRowsEngineFactory::~xmlFilterRowsEngineFactory()
{
    ///@TODO: delete trees
}
