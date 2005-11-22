/**@file xmlCreateColumnsEngineFactory.cxx

@brief implementation of class xmlCreateColumnsEngineFactory

$Header$
*/

#include "xmlCreateColumnsEngineFactory.h"
#include "src/ImActivityNodes/CreateColumnsEngineNode.h"
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

xmlCreateColumnsEngineFactory::xmlCreateColumnsEngineFactory(XTExprsnParser& parser) : xmlFactoryBase(parser)
{
}

IImActivityNode* xmlCreateColumnsEngineFactory::operator()(const DOMElement* xmlActivityNode)
{
    // Retrieve name and node id
    DOMElement* displayInfo = xmlBase::Dom::findFirstChildByName(xmlActivityNode, "DisplayInfo");
    std::string sType       = "CreateColumnsEngineNode";
    std::string sName       = xmlBase::Dom::getAttribute(displayInfo, "labelText");
    std::string sId         = xmlBase::Dom::getAttribute(xmlActivityNode, "id");

    // Create the node
    CreateColumnsEngineNode* node = new CreateColumnsEngineNode(sType, sName, sId);

    // Need to find the list of variables...
    StringList columnNames;
    StringList colExpressions;
    std::vector<StringList> parsedColExps;
    CreateColumnsEngineNode::XprsnNodeMap xprsnNodeMap;

    DOMEvector xmlColumnsVec = getXTSubPropVec(xmlActivityNode, "newColumns");

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

        // Look to see what type of variable is being created here
        std::string sVarType = xmlBase::Dom::getAttribute(xmlColVarVec[0], "value");
        if (sVarType == "categorical") continue;

        columnNames.push_back(sVarName);

        std::string sExpression = "";

        try
        {
            DOMElement* xmlComplex = xmlBase::Dom::findFirstChildByName(xmlColVarVec[1], "Complex");
            if (xmlComplex == 0) throw XTENexception("xmlCreateColumnsFactory finds zero pointer to xmlComplex");
            sExpression = xmlBase::Dom::getTextContent(xmlComplex);
        }
        //catch(xmlBase::WrongNodeType& e)
        catch(...)
        {
            sExpression = xmlBase::Dom::getAttribute(xmlColVarVec[1], "value");
        }
            
        colExpressions.push_back(sExpression);

        // Parse
        IXTExprsnNode* xprsn = XprsnParser().parseExpression(sExpression);

        // Get the tuple column value pointer
        XTcolumnVal<double>* xtColumnVal = 0;
        XTcolumnVal<double>::XTtupleMap::iterator dataIter = XprsnParser().getXtTupleVars().find(sVarName);
            
        if (dataIter != XprsnParser().getXtTupleVars().end()) xtColumnVal = dataIter->second;
        else
        {
            xtColumnVal = new XTcolumnVal<double>(sVarName);
            XprsnParser().getXtTupleVars()[sVarName] = xtColumnVal;
        }

        xprsnNodeMap[xtColumnVal] = xprsn;
    }

    // Set the expressions
    if (!columnNames.empty())
    {
        node->setColumnNames(columnNames);
        node->setXprsnNodeMap(xprsnNodeMap);
    }

    return node;
}

xmlCreateColumnsEngineFactory::~xmlCreateColumnsEngineFactory()
{
    ///@TODO: delete trees
}
