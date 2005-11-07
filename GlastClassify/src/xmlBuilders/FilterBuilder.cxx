/** @file DecisionTreeBuilder.cxx
 *    @brief implementation of classification::Tree; declaration and implementation or its private helper classification::DecisionTreeBuilder::Node
 *
 *    $Header$
 */

#include "DecisionTreeBuilder.h"
#include "xmlBase/XmlParser.h"
#include "facilities/Util.h"
#include <xercesc/dom/DOMElement.hpp>
#include "xmlBase/Dom.h"

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include <set> 
using std::cout;
using std::cerr;
using std::endl;

namespace {

    // Helper function creates a std::string from a DOMString
    // Usage:
    // a = getNextWord(sList, iDelimPos=-1);
    // b = getNextWord(sList, iDelimPos);

    std::string getNextWord(std::string &sList, int &iEnd)
    {
        // This allows space delimited data as well as comma, semicolon, etc.
        static char  cDELIM = ' '; 
        const char *cList;
        int iStart = iEnd + 1; // Start is previous end plus one.
        if(iStart > (int)sList.length())
        {
            // Serious error when you don't initialize iDelimPos=-1
            // before calling this function!
            cerr << "Error with getNextWord!" << endl;
            return sList;
        }
        cList = sList.c_str();
        while((cList[iStart] == cDELIM) && cList[iStart] != 0)
        { 
            iStart++; 
        }
        iEnd = sList.find(cDELIM, iStart);
        if(iEnd > (int)sList.length())
        {
            cerr << "Error with getNextWord" << endl;
        }
        return sList.substr(iStart, iEnd - iStart);
    }

    double getNextDouble(std::string &sList, int &iEnd)
    {
        std::string sIValue = getNextWord(sList, iEnd);
        return atof(sIValue.c_str());
    }


    // output errors to console. (note: m_ prefex for data members )
    // only for statistics to show in debug output
    int  nesting_level, max_depth, node_count;  
    double min_prob, max_prob;
    
    std::string indent()
    {
        std::string ret("  ");; 
        for( int i =0; i<nesting_level;++i) ret += "   ";
        return ret;
    }
} // anonomous namespace

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
//namespace classifier
//{

XERCES_CPP_NAMESPACE_USE

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/*
* The idea for classification::Tree is to encapuslate the xml so that the
*   user only
*   has to provide a filename or string and we will give them access to the
*   linked list of the data which will have the ability to do neural 
*   networking algorithms.
*/


DecisionTreeBuilder::DecisionTreeBuilder(const DOMDocument* document, 
                                         DecisionTreeBuilder::ILookUpIndex& searcher, 
                                         std::ostream& log, 
                                         int iVerbosity )
  : m_domDocument(document),
    m_searcher(searcher),
    m_log(log),
    m_outputLevel(iVerbosity)

  //:m_fields(0) //unless overriden use old method for access to data values.
{
    
    if(m_domDocument == 0)
    {
        // Error checking usually for a missing file.
        // Or sometimes due to a bad XML document.
        // Remember <P><BR></P> is only valid
        // in poorly formed html, not xml.
        // When we get a schema for UserLibrary, we'll
        // be able to use that as a validator also.
        throw Exception( "Error: invalid input file ");
    }

    // Clear the variable name list
    m_varNames.clear();

    // Tree(-1) gives us no stdout.
    // Tree(0) and Tree() gives us errors only
    // Tree(1) gives us errors and verbose output.
    nesting_level = 0;
    max_depth     = 0;
    node_count    = 0;
    min_prob      = 1e99; 
    max_prob      =-1e99;
  
    return;
}


DecisionTreeBuilder::~DecisionTreeBuilder()
{
    return;
}
  
DecisionTree* DecisionTreeBuilder::buildTree(const std::string &treeName)
{
    // Default return value
    DecisionTree* decisionTree = 0;

    m_varNames.clear();
    
    DOMElement* domRoot = m_domDocument->getDocumentElement();

    // This assumes that the node under the root node is Activity Node.
    // In IMW files, this is not the case.
    // IMML/Worksheet/ActivityNodeList
    std::vector<std::string> straListPath;
    straListPath.push_back("Worksheet");
    straListPath.push_back("ActivityNodeList");
    DOMElement* xmlList = findXPath(domRoot, straListPath);
    // If they are using the user library instead of IMW file, this will work.
    if(xmlList == 0) xmlList = domRoot;

    // We now need to ensure that we're getting 
    // ActivityNode[@engineClass=='com.insightful.miner.PredictEngineNode']
    std::vector<DOMElement *> xmlActivityNodes;
    xmlBase::Dom::getChildrenByTagName(xmlList, "ActivityNode", xmlActivityNodes);

    unsigned int iNode=0;
    if(m_outputLevel>1)
    {
        m_log << "Prediction Nodes:\n" 
              << std::setw(10) << "id"    << std::setw(30) << "label" 
              << std::setw(10) << "nodes" << std::setw(10) << "max depth" 
              << std::setw(12) << "min prob" 
              << std::setw(12) << "max prob" << std::endl;
    }

    // Create path to TreeList objects in the node
    std::vector<std::string> path2TreeList;
    path2TreeList.push_back("ModelProperties");
    path2TreeList.push_back("IMML");
    path2TreeList.push_back("TreeList");
    //straModelPath2.push_back("TreeModel");
    
    for(iNode = 0; iNode < xmlActivityNodes.size(); iNode++)
    {
        DOMElement* xmlActivityNode = xmlActivityNodes[iNode];

        std::string sType = xmlBase::Dom::getAttribute(xmlActivityNode, "engineClass");  

        //Look for "PredictEngineNode"
        if (sType == "com.insightful.miner.PredictEngineNode")
        { 
            // Retrieve the name of this tree
            DOMElement* displayInfo = xmlBase::Dom::findFirstChildByName(xmlActivityNode, "DisplayInfo");
            std::string node_name   = xmlBase::Dom::getTagName(displayInfo);
            std::string label_name  = xmlBase::Dom::getAttribute(displayInfo, "labelText");

            // if not the one we want then skip
            if (label_name != treeName)
            {
                continue;
            }

            // Recover the output variable name
            std::string m_outVarName = getCTOutputName(xmlActivityNode);

            // Find the "TreeList" element
            DOMElement* xmlTreeList = findXPath(xmlActivityNode, path2TreeList);

            // Build the independent variable map
            tupleVarIndexMap tupleVarMap = buildVarIndexMap(xmlTreeList);

            // Retrieve all instances of "TreeModel" in this "TreeList"
            std::vector<DOMElement *> xmlTreeModelVec;
            xmlBase::Dom::getChildrenByTagName(xmlTreeList, "TreeModel", xmlTreeModelVec);

            // This canna happen so must be worthy of an exception
            if (xmlTreeModelVec.empty())
            {
                throw Exception("ActivityNode/ModelProperties/IMML/TreeList/TreeModel not found.");
            }

            // Create a new DecisionTree object
            decisionTree = new DecisionTree(label_name);

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

            // Presumably, we are done
            break;
        }
    }

    return decisionTree;
}


  // Inputs a node (usually the root node)
  // and a list such as:
  //	char *straXMLPath[] = {"ActivityNode", "ModelProperties", "IMML", 
  //                           "TreeList", "TreeModel", 0};
  // And returns the node at the end, in this case: TreeModel.
  // The node given should be:
  // <xmlParent><ActivityNode><ModelProperties>
  //		<IMML><TreeList><TreeModel></TreeModel></TreeList></IMML>
  // </ModelProperties></ActivityNode></xmlParent>
  // When calling this function be sure to check Node==0 
  // or Node.isNull()==true;
  // Replace above function entirely because 
  //  - unnecessary calls to native Xerces routines
  //  - vector of strings is much nicer to deal with than array
  //  - all this requires is a simple loop, not recursion
DOMElement* DecisionTreeBuilder::findXPath(DOMElement* xmlParent, const std::vector<std::string>& nodeNames)
{
    unsigned nNames = nodeNames.size();
    if (nNames == 0) return xmlParent;

    DOMElement* parent = xmlParent;
    DOMElement* child;
    for (unsigned iName = 0; iName < nNames; iName++) 
    {
        child = xmlBase::Dom::findFirstChildByName(parent, nodeNames[iName]);
        if (child == 0) return child;  // path ended prematurely
        parent = child;
    }
    
    return child;
}
  // This function takes a TreeModel Root node and parses its children
  // It currently only works for UserLibrary.xml and may not be compatible
  // with multiple outputs. It will in the future. All it does is call 
  // parseFields and parseNode. To handle multiple outputs, I would assume 
  // that one would call parseNode for all TreeModel's Node children. However,
  // this would require multiple Node trees, one for each.
  // This could be handled by having the root have a child for each node, but 
  // I am going to wait until I get a better view of what is going on.

void DecisionTreeBuilder::parseTree(DOMElement* xmlTreeModel, DecisionTree* decisionTree, std::map<std::string, int>& varIndexMap)
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
    parseNode(xmlFirstNode, decisionTree, varIndexMap);
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

void DecisionTreeBuilder::parseNode(DOMElement* xmlElement, DecisionTree* decisionTree, std::map<std::string, int>& varIndexMap)
{
    // Set the depth and node count
    ++nesting_level; 
    ++node_count;
    max_depth = max_depth> nesting_level? max_depth : nesting_level;

    // determine the node ID
    std::string sID = xmlBase::Dom::getAttribute(xmlElement, "id");
    if(m_outputLevel>2) m_log <<  indent() <<  "Node " << sID << ": ";
    int nodeId = atoi(sID.c_str());

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
            double weight = 0.;

            // Get the weight for this node... will depend on whether a classification tree
            // or a regression tree. 
            std::string sList = xmlBase::Dom::getAttribute(xmlElement, "yprob");
            std::string score = xmlBase::Dom::getAttribute(xmlElement, "score");
            if( sList.empty()) 
            {
                // must be a regression tree prediction node: use the probability list, just one element for the value
                weight = facilities::Util::stringToDouble(score);
            }
            else
            {
                if( m_outputLevel>2 ) m_log << "yprob: '" << sList << "';" << endl;
                int iDelimPos = -1;
                //while((iDelimPos > -1) ^ (m_yprob.size() == 0))
                //{ 
                //    // iDelimPos == -1 means that we've gone around the list. 
                //    // m_yprob.size() == 0 means that there's nothing in the list.
                //    // Beware, if you get an endless loop, here it lies.
                //    m_yprob.push_back(getNextDouble(sList, iDelimPos));
                //}
                // Not completely sure what to do here... take the first value and proceed
                weight = getNextDouble(sList, iDelimPos);
            }

            // Create the node
            decisionTree->addNode(nodeId, -1, weight);

            min_prob = std::min(min_prob, weight);
            max_prob = std::max(max_prob, weight);
        }
    }
    // Otherwise we are at a branch point
    else
    {
        std::string varname   = xmlBase::Dom::getAttribute(xmlPredicate, "field");
        int         index     = varIndexMap[varname];
        double      value     = xmlBase::Dom::getDoubleAttribute(xmlPredicate, "value");
        std::string sOperator = xmlBase::Dom::getAttribute(xmlPredicate, "operator");

        // Add a new node
        decisionTree->addNode(nodeId, index, value);

        // how many nodes?
        if (xmlNodeVec.size() != 2)
        {
            int j = xmlNodeVec.size();
            int i = 0;
        }

        // Order of children depends on the operator
        if (sOperator == "lessThan")
        {
            parseNode(xmlNodeVec[0], decisionTree, varIndexMap);
            parseNode(xmlNodeVec[1], decisionTree, varIndexMap);
        }
        else
        {
            parseNode(xmlNodeVec[1], decisionTree, varIndexMap);
            parseNode(xmlNodeVec[0], decisionTree, varIndexMap);
        }
    }

    --nesting_level;

    return;
}

std::string DecisionTreeBuilder::getCTOutputName(DOMElement* xmlActivityNode)
{
    std::string output;

    // Create path to Argument output list in the node
    std::vector<std::string> path2Arguments;
    path2Arguments.push_back("ArgumentList");
    path2Arguments.push_back("XTProps");

    // Search for the output argument list, this to find the name of the output of this CT
    DOMElement* argumentList = findXPath(xmlActivityNode, path2Arguments);

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

DecisionTreeBuilder::tupleVarIndexMap DecisionTreeBuilder::buildVarIndexMap(DOMElement* xmlTreeList)
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


//}