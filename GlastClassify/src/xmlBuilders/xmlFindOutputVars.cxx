/**@file xmlFindOutputVars.cxx

@brief implementation of class xmlFindOutputVars

$Header$
*/

#include "xmlFindOutputVars.h"
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

xmlFindOutputVars::xmlFindOutputVars(XTExprsnParser& parser) : xmlFactoryBase(parser)
{
}

int xmlFindOutputVars::operator()(const DOMElement* xmlActivityNode)
{
    int numVars = 0;

    // Retrieve name and node id
    DOMElement* displayInfo = xmlBase::Dom::findFirstChildByName(xmlActivityNode, "DisplayInfo");
    std::string sType       = "CreateColumnsEngineNode";
    std::string sName       = xmlBase::Dom::getAttribute(displayInfo, "labelText");
    std::string sId         = xmlBase::Dom::getAttribute(xmlActivityNode, "id");

    // Need to find the list of variables...
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

        // Get the tuple column value pointer
        XTcolumnValBase*     xtColumnVal = 0;
        XTtupleMap::iterator dataIter = XprsnParser().getXtTupleVars().find(sVarName);
            
        if (dataIter == XprsnParser().getXtTupleVars().end())
        {
            numVars++;

            if (sVarType == "continuous") xtColumnVal = new XTcolumnVal<REALNUM>(sVarName);
            else                          xtColumnVal = new XTcolumnVal<std::string>(sVarName,"categorical");

            XprsnParser().getXtTupleVars()[sVarName] = xtColumnVal;
        }
    }

    return numVars;
}

xmlFindOutputVars::~xmlFindOutputVars()
{
    ///@TODO: delete trees
}
