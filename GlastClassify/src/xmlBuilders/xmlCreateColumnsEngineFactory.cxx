/**@file xmlCreateColumnsEngineFactory.cxx

@brief implementation of class xmlCreateColumnsEngineFactory

$Header$
*/

#include "xmlCreateColumnsEngineFactory.h"
#include "../ImActivityNodes/CreateColumnsEngineNode.h"
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

xmlCreateColumnsEngineFactory::xmlCreateColumnsEngineFactory(std::ostream& log, int iVerbosity )
                        : xmlFactoryBase(log,iVerbosity), 
                          m_log(log), m_outputLevel(iVerbosity)
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

    DOMEvector xmlColumnsVec = getXTSubPropVec(xmlActivityNode, "newColumns");

    for(DOMEvector::iterator colPropIter = xmlColumnsVec.begin(); colPropIter != xmlColumnsVec.end(); colPropIter++)
    {
        DOMElement* xmlColVar = *colPropIter;

        DOMEvector xmlColVarVec;
        xmlBase::Dom::getChildrenByTagName(xmlColVar, "Property", xmlColVarVec);

        std::string sVarName    = xmlBase::Dom::getAttribute(xmlColVar, "name");
        std::string sExpression = "";

        columnNames.push_back(sVarName);

        try
        {
            DOMElement* xmlComplex = xmlBase::Dom::findFirstChildByName(xmlColVarVec[1], "Complex");
            sExpression = xmlBase::Dom::getTextContent(xmlComplex);
        }
        //catch(xmlBase::WrongNodeType& e)
        catch(...)
        {
            sExpression = xmlBase::Dom::getAttribute(xmlColVarVec[1], "value");
        }

        // Trim the blank spaces
        sExpression = trimBlanks(sExpression);
            
        colExpressions.push_back(sExpression);

        // Parse
        StringList parsedExpression;
        parseExpression(parsedExpression, sExpression);

        parsedColExps.push_back(parsedExpression);
    }

    // Set the expressions
    if (!columnNames.empty())
    {
        node->setColumnNames(columnNames);
        node->setColumnExpressions(colExpressions);
        node->setParsedColExps(parsedColExps);
    }

    return node;
}

xmlCreateColumnsEngineFactory::~xmlCreateColumnsEngineFactory()
{
    ///@TODO: delete trees
}
