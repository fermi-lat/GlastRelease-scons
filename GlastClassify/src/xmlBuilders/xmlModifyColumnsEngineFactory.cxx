/**@file xmlModifyColumnsEngineFactory.cxx

@brief implementation of class xmlModifyColumnsEngineFactory

$Header$
*/

#include "xmlModifyColumnsEngineFactory.h"
#include "../ImActivityNodes/ModifyColumnsEngineNode.h"
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

xmlModifyColumnsEngineFactory::xmlModifyColumnsEngineFactory(XTExprsnParser& parser) : xmlFactoryBase(parser)
{
}

IImActivityNode* xmlModifyColumnsEngineFactory::operator()(const DOMElement* xmlActivityNode)
{
    // Retrieve name and node id
    DOMElement* displayInfo = xmlBase::Dom::findFirstChildByName(xmlActivityNode, "DisplayInfo");
    std::string sType       = "ModifyColumnsEngineNode";
    std::string sName       = xmlBase::Dom::getAttribute(displayInfo, "labelText");
    std::string sId         = xmlBase::Dom::getAttribute(xmlActivityNode, "id");

    // Create the node
    ModifyColumnsEngineNode* node = new ModifyColumnsEngineNode(sType, sName, sId);

    return node;
}

xmlModifyColumnsEngineFactory::~xmlModifyColumnsEngineFactory()
{
    ///@TODO: delete trees
}
