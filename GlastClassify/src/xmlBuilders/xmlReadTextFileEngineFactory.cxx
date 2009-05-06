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

    // In this engine, we need to make the connection to any input ntuple variables so they can
    // be accessed during the running of the analysis. 
    DOMEvector xmlColumnsVec = getXTSubPropVec(xmlActivityNode, "columns");

    // Loop through the "columns" 
    for(DOMEvector::iterator colPropIter = xmlColumnsVec.begin(); colPropIter != xmlColumnsVec.end(); colPropIter++)
    {
        DOMElement* xmlColVar = *colPropIter;

        DOMEvector xmlColVarVec;
        xmlBase::Dom::getChildrenByTagName(xmlColVar, "Property", xmlColVarVec);

        std::string sVarName = xmlBase::Dom::getAttribute(xmlColVar, "name");

        // Let us check to see if the variable has already been set up (which should not be the case)
        XTcolumnValBase*     xtColumnVal = 0;
        XTtupleMap::iterator dataIter = XprsnParser().getXtTupleVars().find(sVarName);

        // Not found means we create the entry
        if (dataIter == XprsnParser().getXtTupleVars().end())
        {
            // May 6, 2009 - Bill Atwood asks that we not create new variables in this node
            //               in order to not include now obsolete CTB variables

            // For now we assume that all input variables are going to be "continuous" 
            //xtColumnVal = new XTcolumnVal<REALNUM>(sVarName);

            // Add to the list...
            //XprsnParser().getXtTupleVars()[sVarName] = xtColumnVal;
        }
    }

    // Create a flag for determining whether to keep a row in the end of processing
    // @TODO need to change this to a bool value (implement storage maps for bool and categorical vars)
    std::string           sVarName    = "WriteTupleRow";
    XTcolumnValBase*      basePtr     = XprsnParser().getXtTupleVars()[sVarName];
    XTcolumnVal<REALNUM>* xtColumnVal = dynamic_cast<XTcolumnVal<REALNUM>*>(basePtr);

    if (xtColumnVal == 0)
    {
        xtColumnVal = new XTcolumnVal<REALNUM>(sVarName);
        XprsnParser().getXtTupleVars()[sVarName] = xtColumnVal;
    }

    node->setXtColumnVal(xtColumnVal);

    return node;
}

xmlReadTextFileEngineFactory::~xmlReadTextFileEngineFactory()
{
    ///@TODO: delete trees
}
