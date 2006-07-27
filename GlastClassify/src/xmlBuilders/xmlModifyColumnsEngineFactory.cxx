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

    // Need to find the list of variables...
    StringList columnNames;
    ModifyColumnsEngineNode::XprsnNodeMap xprsnNodeMap;

    // Look for "columns" which we are going to modify
    DOMEvector xmlColumnsVec = getXTSubPropVec(xmlActivityNode, "columns");

    for(DOMEvector::iterator colPropIter = xmlColumnsVec.begin(); colPropIter != xmlColumnsVec.end(); colPropIter++)
    {
        DOMElement* xmlColVar = *colPropIter;

        DOMEvector xmlColVarVec;
        xmlBase::Dom::getChildrenByTagName(xmlColVar, "Property", xmlColVarVec);

        std::string sVarName = xmlBase::Dom::getAttribute(xmlColVar, "name");

        // Special code to strip off "Pr(" from the name
        int prPos = sVarName.find("Pr(",0);
        if (prPos > -1)
        {
            sVarName = sVarName.substr(prPos+3, sVarName.length()-prPos-4);
        }

        // Get the information on the type of modification being performed
        std::string sModType = xmlBase::Dom::getAttribute(xmlColVarVec[0], "name");
        std::string sNewName = xmlBase::Dom::getAttribute(xmlColVarVec[0], "value");

        // As of right now (6/30/2006) we can only handle assignment to existing variables...
        if (sModType != "newName") continue;

        // Retrieve the original variable name from the tuple list
        XTtupleMap::iterator dataIter = XprsnParser().getXtTupleVars().find(sVarName);

        // If it doesn't exist then there is a problem (and we should throw an exception...)
        if (dataIter == XprsnParser().getXtTupleVars().end()) continue;

        // Set pointer to the "old" value
        XTcolumnValBase* xtOldValue  = dataIter->second;

        // Look up the "new" value in case it exists
        dataIter = XprsnParser().getXtTupleVars().find(sNewName);
        XTcolumnValBase*     xtColumnVal = 0;
            
        // Exists, get pointer to it
        if (dataIter != XprsnParser().getXtTupleVars().end())
        {
            xtColumnVal = dataIter->second;
        }
        // Doesn't exist, create new one
        else
        {
            std::string sVarType = xtOldValue->getType();

            // New type depends on the old type
            if (sVarType == "continuous") xtColumnVal = new XTcolumnVal<REALNUM>(sVarName);
            else                          xtColumnVal = new XTcolumnVal<std::string>(sVarName,"categorical");

            XprsnParser().getXtTupleVars()[sVarName] = xtColumnVal;
        }

        xprsnNodeMap[xtColumnVal] = xtOldValue;
    }

    // Set the expressions
    if (!xprsnNodeMap.empty())
    {
        node->setXprsnNodeMap(xprsnNodeMap);
    }

    return node;
}

xmlModifyColumnsEngineFactory::~xmlModifyColumnsEngineFactory()
{
    ///@TODO: delete trees
}
