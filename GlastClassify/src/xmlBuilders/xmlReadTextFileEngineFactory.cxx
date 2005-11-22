/**@file xmlReadTextFileEngineFactory.cxx

@brief implementation of class xmlReadTextFileEngineFactory

$Header$
*/

#include "xmlReadTextFileEngineFactory.h"
#include "../ImActivityNodes/ReadTextFileEngineNode.h"
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

xmlReadTextFileEngineFactory::xmlReadTextFileEngineFactory(XTExprsnParser& parser) : xmlFactoryBase(parser)
{
}

IImActivityNode* xmlReadTextFileEngineFactory::operator()(const DOMElement* xmlActivityNode)
{
    // Retrieve name and node id
    DOMElement* displayInfo = xmlBase::Dom::findFirstChildByName(xmlActivityNode, "DisplayInfo");
    std::string sType       = "ReadTextFileEngineNode";
    std::string sName       = xmlBase::Dom::getAttribute(displayInfo, "labelText");
    std::string sId         = xmlBase::Dom::getAttribute(xmlActivityNode, "id");

    // Create the node
    ReadTextFileEngineNode* node = new ReadTextFileEngineNode(sType, sName, sId);

    // Create a flag for determining whether to keep a row in the end of processing
    // @TODO need to change this to a bool value (implement storage maps for bool and categorical vars)
    std::string sVarName = "WriteTupleRow";
    XTcolumnVal<double>* xtColumnVal = new XTcolumnVal<double>(sVarName);
    XprsnParser().getXtTupleVars()[sVarName] = xtColumnVal;

    node->setXtColumnVal(xtColumnVal);

    return node;
}

xmlReadTextFileEngineFactory::~xmlReadTextFileEngineFactory()
{
    ///@TODO: delete trees
}
