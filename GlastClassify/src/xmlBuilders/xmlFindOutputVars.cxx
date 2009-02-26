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
    std::string sName       = xmlBase::Dom::getAttribute(displayInfo, "labelText");
    std::string sId         = xmlBase::Dom::getAttribute(xmlActivityNode, "id");

    std::string sType       = xmlBase::Dom::getAttribute(xmlActivityNode, "engineClass");

    int idx = sType.find("CreateColumns");
    if (idx > -1)
    {
        numVars = numCreateColumnsVars(xmlActivityNode);
    }
    else if ((idx = sType.find("PredictEngine")) > -1)
    {
        numVars = numPredictEngineVars(xmlActivityNode);
    }

    return numVars;
}

int xmlFindOutputVars::numCreateColumnsVars(const DOMElement* xmlActivityNode)
{
    // Searches input document for variables created by "CreateColumns" nodes
    // Stores the found variables in the "local ntuple" list for access by 
    // other nodes during processing 
    int numVars = 0;

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

int xmlFindOutputVars::numPredictEngineVars(const DOMElement* xmlActivityNode)
{
    // Searches the input document for output variables created by the "Predict" nodes
    int numVars = 0;

    // Create path to Argument output list in the node
    std::vector<std::string> path2Arguments;
    path2Arguments.push_back("ArgumentList");
    path2Arguments.push_back("XTProps");

    // Search for the output argument list, this to find the name of the output of this CT
    const DOMElement* argumentList = findXPath(xmlActivityNode, path2Arguments);

    // Obtain the list of properties
    std::vector<DOMElement*> xmlArgListVec;
    xmlBase::Dom::getChildrenByTagName(argumentList, "Property", xmlArgListVec);
        
    // Search through the list of properties for the "NewColumns" property 
    // (which is probably the first property?)
    for(std::vector<DOMElement*>::iterator xmlArgListVecItr = xmlArgListVec.begin();
        xmlArgListVecItr != xmlArgListVec.end(); xmlArgListVecItr++)
    {
        DOMElement* xmlArgList = *xmlArgListVecItr;
           
        std::string argListTag  = xmlBase::Dom::getTagName(xmlArgList);
        std::string argName     = xmlBase::Dom::getAttribute(xmlArgList, "name");

        // "NewColumns"?
        if (argName == "newColumns")
        {
            // Get this list of properties here
            std::vector<DOMElement*> xmlPropListVec;
            xmlBase::Dom::getChildrenByTagName(xmlArgList, "Property", xmlPropListVec);
    
            for(std::vector<DOMElement*>::iterator xmlPropListVecItr = xmlPropListVec.begin();
                xmlPropListVecItr != xmlPropListVec.end(); xmlPropListVecItr++)
            {
                DOMElement* xmlProperty = *xmlPropListVecItr;
           
                std::string propertyTag  = xmlBase::Dom::getTagName(xmlProperty);
                std::string propertyName = xmlBase::Dom::getAttribute(xmlProperty, "name");

                // "NewColumns"?
                if (propertyName == "specifiedCategory")
                {
                    std::string sName = xmlBase::Dom::getAttribute(xmlProperty, "value");
    
                    // But make sure it doesn't already exist!
                    XTtupleMap::iterator dataIter = XprsnParser().getXtTupleVars().find(sName);
  
                    // Doesn't exist, add it
                    if (dataIter == XprsnParser().getXtTupleVars().end())
                    {
                        XTcolumnValBase* xtColumnVal = new XTcolumnVal<REALNUM>(sName);

                        XprsnParser().getXtTupleVars()[sName] = xtColumnVal;

                        numVars++;
                    }
                }
            }

            break;
        }
    }

    return numVars;
}

xmlFindOutputVars::~xmlFindOutputVars()
{
    ///@TODO: delete trees
}
