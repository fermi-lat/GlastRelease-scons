/**@file xmlPredictEngineFactory.cxx

@brief implementation of class xmlPredictEngineFactory

$Header$
*/

#include "xmlPredictEngineFactory.h"
#include "../ImActivityNodes/PredictEngineNode.h"
#include "classifier/DecisionTree.h"
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
    /** @class Exception 
        @brief hold a string
    */
    class Exception : public std::exception
    {
    public: 
        Exception(std::string error):m_what(error){}
        ~Exception() throw() {;}
        virtual const char *what( ) const  throw() { return m_what.c_str();} 
        std::string m_what;
    };

    // output errors to console. (note: m_ prefex for data members )
    // only for statistics to show in debug output
    int  nesting_level, max_depth, node_count;  
    double min_prob, max_prob;
} // anonomous namespace

xmlPredictEngineFactory::xmlPredictEngineFactory(std::ostream& log, int iVerbosity )
                        : xmlFactoryBase(log,iVerbosity), 
                          m_yProbIndex(0), m_numProbVals(0), m_log(log), m_outputLevel(iVerbosity)
{
}

IImActivityNode* xmlPredictEngineFactory::operator()(const DOMElement* xmlActivityNode)
{
    // Retrieve name and node id
    DOMElement* displayInfo = xmlBase::Dom::findFirstChildByName(xmlActivityNode, "DisplayInfo");
    std::string sType       = "PredictEngineNode";
    std::string sName       = xmlBase::Dom::getAttribute(displayInfo, "labelText");
    std::string sId         = xmlBase::Dom::getAttribute(xmlActivityNode, "id");

    // Create the node
    PredictEngineNode* node = new PredictEngineNode(sType, sName, sId);

    DecisionTree* tree = parseForest(xmlActivityNode);

    node->setDecisionTree(tree);
    node->setInputVar(m_varNames);
    node->addOutputVar(m_outVarName);

    return node;
}

xmlPredictEngineFactory::~xmlPredictEngineFactory()
{
    ///@TODO: delete trees
}

DecisionTree* xmlPredictEngineFactory::parseForest(const DOMElement* xmlActivityNode)
{
    // Recover the output variable name
    m_specCatName = getCTOutputName(xmlActivityNode); 
    m_outVarName  = "Pr(" + m_specCatName + ")";

    // Use this to decode the order of probabilities we'll encounter
    std::vector<DOMElement *> xmlColInfoVec;
    xmlBase::Dom::getDescendantsByTagName(xmlActivityNode, "ColumnInfo", xmlColInfoVec);

    m_yProbIndex  = 0;
    m_numProbVals = 0;

    // Loop through Info nodes to find output tags
    for(std::vector<DOMElement*>::iterator colInfoIter = xmlColInfoVec.begin(); 
        colInfoIter != xmlColInfoVec.end(); colInfoIter++)
    {
        DOMElement* xmlColInfo = *colInfoIter;

        std::string sName = xmlBase::Dom::getAttribute(xmlColInfo, "name");
        std::string sRole = xmlBase::Dom::getAttribute(xmlColInfo, "role");

        if (sRole == "dependent")
        {
            DOMElement* xmlLevel = xmlBase::Dom::findFirstChildByName(xmlColInfo, "Level");

            while(xmlLevel != 0)
            {
                std::string sLevelName = xmlBase::Dom::getAttribute(xmlLevel, "value");
                
                if (sLevelName == m_specCatName) m_yProbIndex = m_numProbVals;
                
                xmlLevel = xmlBase::Dom::getSiblingElement(xmlLevel);
                m_numProbVals++;
            }
            break;
        }
    }

    // Build the independent variable map
    //tupleVarIndexMap tupleVarMap = buildVarIndexMap(xmlTreeList);
    tupleVarIndexMap tupleVarMap;

    tupleVarMap.clear();
    m_varNames.clear();

    // Retrieve all instances of "TreeModel" in this "TreeList"
    std::vector<DOMElement *> xmlTreeModelVec;
    xmlBase::Dom::getDescendantsByTagName(xmlActivityNode, "TreeModel", xmlTreeModelVec);

    // This canna happen so must be worthy of an exception
    if (xmlTreeModelVec.empty())
    {
        throw Exception("ActivityNode/ModelProperties/IMML/TreeList/TreeModel not found.");
    }

    // Retrieve the name of this tree again...
    DOMElement* displayInfo = xmlBase::Dom::findFirstChildByName(xmlActivityNode, "DisplayInfo");
    std::string label_name  = xmlBase::Dom::getAttribute(displayInfo, "labelText");

    // Create a new DecisionTree object
    DecisionTree* decisionTree = new DecisionTree(label_name);

    // Create a weight for each tree (based on the number of trees)
    double treeWeight = 1. / xmlTreeModelVec.size();

    // Loop over TreeModels in the TreeList
    for(std::vector<DOMElement*>::iterator xmlTreeModelItr = xmlTreeModelVec.begin();
        xmlTreeModelItr != xmlTreeModelVec.end(); xmlTreeModelItr++)
    {
        DOMElement* xmlTreeModel = *xmlTreeModelItr;

        if( m_outputLevel>0)
        {
            std::string sId = xmlBase::Dom::getAttribute(xmlActivityNode, "id");
            m_log << std::setw(10) << sId << std::setw(30) << label_name ;
        }

        // Create top node for DecisionTree:
        decisionTree->addNode(0, -10, treeWeight);

        // Ok, now parse the decision tree....
        parseTree(xmlTreeModel, decisionTree, tupleVarMap);

        if (m_outputLevel>1) 
        {
            m_log << std::setw(10) << node_count <<std::setw(10) << max_depth 
                  << std::setw(12) << std::setprecision(3)<< min_prob
                  << std::setw(12) << std::setprecision(3)<< max_prob << std::endl;
        }

        node_count = 0;
        max_depth  = 0;
        min_prob   = 10e99; 
        max_prob   = 0;
    }

    int numTreeModels = xmlTreeModelVec.size();

    if (m_outputLevel > 1)
    {
        m_log << "Number of trees: " << numTreeModels << std::endl;
    }

    return decisionTree;
}

  // This function takes a TreeModel Root node and parses its children
  // It currently only works for UserLibrary.xml and may not be compatible
  // with multiple outputs. It will in the future. All it does is call 
  // parseFields and parseNode. To handle multiple outputs, I would assume 
  // that one would call parseNode for all TreeModel's Node children. However,
  // this would require multiple Node trees, one for each.
  // This could be handled by having the root have a child for each node, but 
  // I am going to wait until I get a better view of what is going on.

void xmlPredictEngineFactory::parseTree(DOMElement* xmlTreeModel, DecisionTree* decisionTree, std::map<std::string, int>& varIndexMap)
{
    // xmlMiningSchema == <MiningSchema>
    //				<MiningField name=""/>
    //				<MiningField name=""/>
    //				...
    //			  </MiningSchema>
    if( !xmlTreeModel->hasChildNodes() )   
    {
        throw Exception(" <TreeModel> is not correct.");	
    }
    // xmlTreeModel == <TreeModel>
    // xmlFirstNode ==  <Node>
    //	     <SimplePredicate field="" operator="" value=""/>
    //	     <Node>
    //		<SimplePredicate field="" operator="" value=""/>
    //		<True />
    //	     </Node>
    //	     <Node>...</Node>
    //	   </Node>
    //	</TreeModel>
    DOMElement* xmlFirstNode = xmlBase::Dom::getFirstChildElement(xmlTreeModel);

    // determine the node ID (done this way to comply with Toby's DecisionTree class...
    std::string sID = xmlBase::Dom::getAttribute(xmlFirstNode, "id");
    if(m_outputLevel>2) m_log <<  indent(nesting_level) <<  "Node " << sID << ": ";
    int nodeId = atoi(sID.c_str());

    parseNode(xmlFirstNode, nodeId, decisionTree, varIndexMap);
}

  // This is a recursive function that parses xmlElement into oNode
  // as well as parsing xmlElement's children into oNode's children.
  /* example log with verbose:
     Node: id: '15';
     Child: SimplePredicate
     Field: 'TKR.1.E1stChisq', 6;
     Value: '2.87719';
     Operator == greaterOrEqual
  */

void xmlPredictEngineFactory::parseNode(DOMElement* xmlElement, int nodeId, DecisionTree* decisionTree, std::map<std::string, int>& varIndexMap)
{
    // Set the depth and node count
    ++nesting_level; 
    ++node_count;
    max_depth = max_depth> nesting_level? max_depth : nesting_level;
/*
    // determine the node ID
    std::string sID = xmlBase::Dom::getAttribute(xmlElement, "id");
    if(m_outputLevel>2) m_log <<  indent() <<  "Node " << sID << ": ";
    int nodeId = atoi(sID.c_str());
*/
    // The first child of the node will be the "predicate" describing the node
    DOMElement* xmlPredicate = xmlBase::Dom::getFirstChildElement(xmlElement);
    std::string sNode = xmlBase::Dom::getTagName(xmlPredicate);

    // Get the list of "Nodes" below this one
    std::vector<DOMElement *> xmlNodeVec;
    xmlBase::Dom::getChildrenByTagName(xmlElement, "Node", xmlNodeVec);

    // If no nodes then we are at a leaf
    if (xmlNodeVec.empty())
    {
        if (sNode == "True")
        {
            // Default value for weight
            double weight[4] = {0.,0.,0.,0.}; // plenty of length for now...

            // Get the weight for this node... will depend on whether a classification tree
            // or a regression tree. 
            std::string sList = xmlBase::Dom::getAttribute(xmlElement, "yprob");
            std::string score = xmlBase::Dom::getAttribute(xmlElement, "score");
            if( sList.empty()) 
            {
                // must be a regression tree prediction node: use the probability list, just one element for the value
                weight[0] = facilities::Util::stringToDouble(score);
            }
            else
            {
                if( m_outputLevel>2 ) m_log << "yprob: '" << sList << "';" << std::endl;
                int iDelimPos = -1;
                //while((iDelimPos > -1) ^ (m_yprob.size() == 0))
                //{ 
                //    // iDelimPos == -1 means that we've gone around the list. 
                //    // m_yprob.size() == 0 means that there's nothing in the list.
                //    // Beware, if you get an endless loop, here it lies.
                //    m_yprob.push_back(getNextDouble(sList, iDelimPos));
                //}
                // Not completely sure what to do here... take the first value and proceed
                int idx = 0;

                weight[idx++] = getNextDouble(sList, iDelimPos);

                // Try the second value if it exists...
                while(iDelimPos > -1) weight[idx++] = getNextDouble(sList, iDelimPos);
            }

            // Create the node
            decisionTree->addNode(nodeId, -1, weight[m_yProbIndex]);

            min_prob = std::min(min_prob, weight[m_yProbIndex]);
            max_prob = std::max(max_prob, weight[m_yProbIndex]);
        }
    }
    // Otherwise we are at a branch point
    else
    {
        std::string varname   = xmlBase::Dom::getAttribute(xmlPredicate, "field");
        double      value     = xmlBase::Dom::getDoubleAttribute(xmlPredicate, "value");
        std::string sOperator = xmlBase::Dom::getAttribute(xmlPredicate, "operator");

        // Check to see if the variable is in our map already
        if (varIndexMap.find(varname) == varIndexMap.end())
        {
            // New variable so add to the list of variables and add to map
            m_varNames.push_back(varname);
            varIndexMap[varname] = varIndexMap.size();
        }

        // Index of variable for the node
        int index = varIndexMap[varname];

        // Add a new node
        decisionTree->addNode(nodeId, index, value);

        // how many nodes?
        //if (xmlNodeVec.size() != 2)
        //{
        //    int j = xmlNodeVec.size();
        //    int i = 0;
        //}

        // Retrieve the node ID's
        std::string sNodeIdLeft = xmlBase::Dom::getAttribute(xmlNodeVec[0], "id");
        int nodeIdLeft = atoi(sNodeIdLeft.c_str());
        std::string sNodeIdRight = xmlBase::Dom::getAttribute(xmlNodeVec[1], "id");
        int nodeIdRight = atoi(sNodeIdRight.c_str());

        // Forget all the above... scheme is:
        nodeIdLeft  = 2 * nodeId;
        nodeIdRight = nodeIdLeft + 1;

        // Order of children depends on the operator
        if (sOperator == "lessThan")
        {
            parseNode(xmlNodeVec[0], nodeIdLeft,  decisionTree, varIndexMap);
            parseNode(xmlNodeVec[1], nodeIdRight, decisionTree, varIndexMap);
        }
        else
        {
            parseNode(xmlNodeVec[1], nodeIdLeft,  decisionTree, varIndexMap);
            parseNode(xmlNodeVec[0], nodeIdRight, decisionTree, varIndexMap);
        }
    }

    --nesting_level;

    return;
}

std::string xmlPredictEngineFactory::getCTOutputName(const DOMElement* xmlActivityNode)
{
    std::string output;

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
                    output = xmlBase::Dom::getAttribute(xmlProperty, "value");
                    break;
                }
            }

            break;
        }
    }

    return output;
}

xmlPredictEngineFactory::tupleVarIndexMap xmlPredictEngineFactory::buildVarIndexMap(DOMElement* xmlTreeList)
{
    // Tuple variable index map
    tupleVarIndexMap tupleVarMap;
    int              index = 0;

    m_varNames.clear();
            
    // Look for the list of variables used by this tree
    DOMElement* xmlVarList = xmlBase::Dom::findFirstChildByName(xmlTreeList, "MiningSchema");

    if (xmlVarList)
    {
        // Can we get a vector of entries for the variables?
        std::vector<DOMElement*> xmlVarListVec;
        xmlBase::Dom::getChildrenByTagName(xmlVarList, "MiningField", xmlVarListVec);
        
        for(std::vector<DOMElement*>::iterator xmlVarListVecItr = xmlVarListVec.begin();
            xmlVarListVecItr != xmlVarListVec.end(); xmlVarListVecItr++)
        {
            DOMElement* xmlVarList = *xmlVarListVecItr;

            std::string varListTag  = xmlBase::Dom::getTagName(xmlVarList);
            std::string varName     = xmlBase::Dom::getAttribute(xmlVarList, "name");
            std::string varListAttr = xmlBase::Dom::getAttribute(xmlVarList, "usageType");

            if (varListAttr == "")
            {
                m_varNames.push_back(varName);

                tupleVarMap[varName] = index++;
            }
        }
    }

    return tupleVarMap;
}
