/**@file xmlSplusGraphEngineFactory.cxx

@brief implementation of class xmlSplusGraphEngineFactory

$Header$
*/

#include "xmlSplusGraphEngineFactory.h"
#include "../ImActivityNodes/TwoDimChartEngineNode.h"
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

xmlSplusGraphEngineFactory::xmlSplusGraphEngineFactory(XTExprsnParser& parser) : xmlFactoryBase(parser)
{
}

IImActivityNode* xmlSplusGraphEngineFactory::operator()(const DOMElement* xmlActivityNode)
{
    // Retrieve name and node id
    DOMElement* displayInfo = xmlBase::Dom::findFirstChildByName(xmlActivityNode, "DisplayInfo");
    std::string sType       = "TwoDimChartEngineNode";
    std::string sName       = xmlBase::Dom::getAttribute(displayInfo, "labelText");
    std::string sId         = xmlBase::Dom::getAttribute(xmlActivityNode, "id");

    // Create the node
    TwoDimChartEngineNode* node = new TwoDimChartEngineNode(sType, sName, sId);

    return node;
}

xmlSplusGraphEngineFactory::~xmlSplusGraphEngineFactory()
{
    ///@TODO: delete trees
}
