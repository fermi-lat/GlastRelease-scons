/**@file xmlFilterColumnsEngineFactory.cxx

@brief implementation of class xmlFilterColumnsEngineFactory

$Header$
*/

#include "xmlFilterColumnsEngineFactory.h"
#include "../ImActivityNodes/FilterColumnsEngineNode.h"
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

xmlFilterColumnsEngineFactory::xmlFilterColumnsEngineFactory(XTExprsnParser& parser) : xmlFactoryBase(parser)
{
}

IImActivityNode* xmlFilterColumnsEngineFactory::operator()(const DOMElement* xmlActivityNode)
{
    // Retrieve name and node id
    DOMElement* displayInfo = xmlBase::Dom::findFirstChildByName(xmlActivityNode, "DisplayInfo");
    std::string sType       = "FilterColumnsEngineNode";
    std::string sName       = xmlBase::Dom::getAttribute(displayInfo, "labelText");
    std::string sId         = xmlBase::Dom::getAttribute(xmlActivityNode, "id");

    // Create the node
    FilterColumnsEngineNode* node = new FilterColumnsEngineNode(sType, sName, sId);

    DOMElement* xmlProperty   = getXTProperty(xmlActivityNode, "excludeColumns");
    DOMElement* xmlExpression = xmlBase::Dom::findFirstChildByName(xmlProperty,  "Property");

    std::string sExpression = xmlBase::Dom::getAttribute(xmlExpression, "name");

    node->setExpression(sExpression);

    return node;
}

xmlFilterColumnsEngineFactory::~xmlFilterColumnsEngineFactory()
{
    ///@TODO: delete trees
}
