/**@file xmlNewPredictEngineFactory.cxx

@brief implementation of class xmlNewPredictEngineFactory

$Header$
*/

#include "xmlNewPredictEngineFactory.h"
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

xmlNewPredictEngineFactory::xmlNewPredictEngineFactory(XTExprsnParser& parser) : xmlFactoryBase(parser),
                                                        m_yProbIndex(0), 
                                                        m_numProbVals(0)
{
    m_catIndex.clear();
}

IImActivityNode* xmlNewPredictEngineFactory::operator()(const DOMElement* xmlActivityNode)
{
    // Retrieve name and node id
    DOMElement* displayInfo = xmlBase::Dom::findFirstChildByName(xmlActivityNode, "DisplayInfo");
    std::string sType       = "NewPredictEngineNode";
    std::string sName       = xmlBase::Dom::getAttribute(displayInfo, "labelText");
    std::string sId         = xmlBase::Dom::getAttribute(xmlActivityNode, "id");

    // Create the node
    newPredictEngineNode* node = new newPredictEngineNode(sType, sName, sId);

    newPredictEngineNode::TreePairVector forest = parseForest(xmlActivityNode);

    // Get the tuple column value pointer
    //XTcolumnVal<double>* xtColumnVal = XprsnParser().getXtTupleVars().addNewDataItem(m_outVarName);
    XTcolumnVal<REALNUM>* xtColumnVal = 0;
    
    XTtupleMap::iterator dataIter = XprsnParser().getXtTupleVars().find(m_outVarName);
            
    if (dataIter != XprsnParser().getXtTupleVars().end())
    {
        XTcolumnValBase* basePtr = dataIter->second;
        
        if (basePtr->getType() == "continuous") xtColumnVal = dynamic_cast<XTcolumnVal<REALNUM>*>(basePtr);
    }
    else
    {
        xtColumnVal = new XTcolumnVal<REALNUM>(m_outVarName);
        XprsnParser().getXtTupleVars()[m_outVarName] = xtColumnVal;
    }

    // This creates a "predict node" output value
    std::string               predClass = "PREDICT.class";
    XTcolumnVal<std::string>* predict   = 0;
    
    dataIter = XprsnParser().getXtTupleVars().find(predClass);
            
    if (dataIter != XprsnParser().getXtTupleVars().end())
    {
        XTcolumnValBase* basePtr = dataIter->second;

        if (basePtr->getType() == "categorical") predict = dynamic_cast<XTcolumnVal<std::string>*>(basePtr);
    }
    else
    {
        predict = new XTcolumnVal<std::string>(predClass, "categorical");
        predict->setDataValue("");
        XprsnParser().getXtTupleVars()[predClass] = predict;
    }

    node->setInputVar(m_varNames);
    node->addOutputVar(m_outVarName);
    node->setXTcolumnVal(xtColumnVal);
    node->setPredictVal(predict);
    node->setTreePairVector(forest);

    return node;
}

xmlNewPredictEngineFactory::~xmlNewPredictEngineFactory()
{
    ///@TODO: delete trees
}

newPredictEngineNode::TreePairVector xmlNewPredictEngineFactory::parseForest(const DOMElement* xmlActivityNode)
{
    // Define the vector to return
    newPredictEngineNode::TreePairVector treePairVec;

    treePairVec.clear();

    // Recover the output variable name
    m_specCatName = getCTOutputName(xmlActivityNode); 
    //m_outVarName  = "Pr(" + m_specCatName + ")";
    m_outVarName  = m_specCatName;

    // Use this to decode the order of probabilities we'll encounter
    std::vector<DOMElement *> xmlColInfoVec;
    xmlBase::Dom::getDescendantsByTagName(xmlActivityNode, "ColumnInfo", xmlColInfoVec);

    m_yProbIndex  = 0;
    m_numProbVals = 0;
    m_catIndex.clear();

    // Loop through Info nodes to find the list of dependent/independent variables
    for(std::vector<DOMElement*>::iterator colInfoIter = xmlColInfoVec.begin(); 
        colInfoIter != xmlColInfoVec.end(); colInfoIter++)
    {
        DOMElement* xmlColInfo = *colInfoIter;

        std::string sName = xmlBase::Dom::getAttribute(xmlColInfo, "name");
        std::string sType = xmlBase::Dom::getAttribute(xmlColInfo, "type");
        std::string sRole = xmlBase::Dom::getAttribute(xmlColInfo, "role");

        // Check to see if this exists in our list already
        XTtupleMap::iterator dataIter = XprsnParser().getXtTupleVars().find(sName);

        // If already there then we can skip
        if (dataIter == XprsnParser().getXtTupleVars().end())
        {
            // Add variable to list
            XTcolumnValBase* basePtr = 0;

            if (sType == "continuous") basePtr = new XTcolumnVal<REALNUM>(sName);
            else
            {
                XTcolumnVal<std::string>* xtColumnVal = new XTcolumnVal<std::string>(sName, "categorical");

                xtColumnVal->setDataValue("");
                basePtr = xtColumnVal;
            }

            XprsnParser().getXtTupleVars()[sName] = basePtr;
        }

        DOMElement* xmlLevel = xmlBase::Dom::findFirstChildByName(xmlColInfo, "Level");

        while(xmlLevel != 0)
        {
            std::string sLevelName = xmlBase::Dom::getAttribute(xmlLevel, "value");
                
            if (sRole == "dependent")
            {
                if (sLevelName == m_specCatName) m_yProbIndex = m_numProbVals;

                m_catIndex[sLevelName] = m_numProbVals;
                
                m_numProbVals++;
            }
        
            // Check to see if this exists in our list already
            XTtupleMap::iterator dataIter = XprsnParser().getXtConstants().find(sLevelName);

            // If already there then we can skip
//            if (dataIter == XprsnParser().getXtTupleVars().end())
            if (dataIter == XprsnParser().getXtConstants().end())
            {
                // Add variable to list
                XTcolumnValBase* basePtr = 0;

                if (sType == "continuous") basePtr = new XTcolumnVal<REALNUM>(sLevelName);
                else
                {
                    XTcolumnVal<std::string>* xtColumnVal = new XTcolumnVal<std::string>(sLevelName, sType);

                    xtColumnVal->setDataValue(sLevelName);
                    basePtr = xtColumnVal;
                }

                XprsnParser().getXtConstants()[sLevelName] = basePtr;
            }
            
            xmlLevel = xmlBase::Dom::getSiblingElement(xmlLevel);
        }
    }
    
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

    // Create a weight for each tree (based on the number of trees)
    double treeWeight = 1. / xmlTreeModelVec.size();

    // Loop over TreeModels in the TreeList
    for(std::vector<DOMElement*>::iterator xmlTreeModelItr = xmlTreeModelVec.begin();
        xmlTreeModelItr != xmlTreeModelVec.end(); xmlTreeModelItr++)
    {
        DOMElement* xmlTreeModel = *xmlTreeModelItr;

        // Ok, now parse the decision tree....
        IXTExprsnNode* treeNode = parseTree(xmlTreeModel);

        newPredictEngineNode::TreePair treePair(treeNode, treeWeight);

        treePairVec.push_back(treePair);

        node_count = 0;
        max_depth  = 0;
        min_prob   = 10e99; 
        max_prob   = 0;
    }

    return treePairVec;
}

  // This function takes a TreeModel Root node and parses its children
  // It currently only works for UserLibrary.xml and may not be compatible
  // with multiple outputs. It will in the future. All it does is call 
  // parseFields and parseNode. To handle multiple outputs, I would assume 
  // that one would call parseNode for all TreeModel's Node children. However,
  // this would require multiple Node trees, one for each.
  // This could be handled by having the root have a child for each node, but 
  // I am going to wait until I get a better view of what is going on.

IXTExprsnNode*  xmlNewPredictEngineFactory::parseTree(DOMElement* xmlTreeModel)
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
    //std::string sID = xmlBase::Dom::getAttribute(xmlFirstNode, "id");
    //int nodeId = atoi(sID.c_str());

    return parseNode(xmlFirstNode);
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

IXTExprsnNode* xmlNewPredictEngineFactory::parseNode(DOMElement* xmlElement)
{
    IXTExprsnNode* node = 0;

    // Set the depth and node count
    ++nesting_level; 
    ++node_count;
    max_depth = max_depth> nesting_level? max_depth : nesting_level;
    
    // determine the node ID
    std::string sNodeId = xmlBase::Dom::getAttribute(xmlElement, "id");
    
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

            // Retrieve all the ancillary info
            std::string sId  = xmlBase::Dom::getAttribute(xmlElement, "id");
            int         id   = 0;
            if (!sId.empty()) id = atoi(sId.c_str());

            std::string sRec = xmlBase::Dom::getAttribute(xmlElement, "recordCount");
            int         rec  = 0;
            if (!sRec.empty()) rec = atoi(sRec.c_str());

            std::string sGrp = xmlBase::Dom::getAttribute(xmlElement, "group");
            int         grp  = 0;
            if (!sGrp.empty()) grp = atoi(sGrp.c_str());

            std::string sDev = xmlBase::Dom::getAttribute(xmlElement, "deviance");
            double      dev  = 0;
            if (!sDev.empty()) dev = atof(sDev.c_str());

            std::string sEnt = xmlBase::Dom::getAttribute(xmlElement, "entropy");
            double      ent  = 0;
            if (!sEnt.empty()) ent = atoi(sEnt.c_str());

            std::string sIni = xmlBase::Dom::getAttribute(xmlElement, "gini");
            double      ini  = 0;
            if (!sIni.empty()) ini = atoi(sIni.c_str());

            std::string sRsk = xmlBase::Dom::getAttribute(xmlElement, "risk");
            int         rsk  = 0;
            if (!sRsk.empty()) rsk = atoi(sRsk.c_str());

            // Create the node
            //double* pValue = new double;

            //*pValue = weight[m_yProbIndex];

            //int         catIdx = m_catIndex[score];
            //double      yProb  = weight[catIdx];
            double      yProb  = weight[m_yProbIndex];
        
            //node = new XTExprsnValue<double>(sNodeId, pValue);
            CTOutPut* ctOutPut = new CTOutPut(id,score,rec,grp,dev,ent,ini,rsk,yProb); 
            node = new XTExprsnValue<CTOutPut>(sNodeId, ctOutPut);

            min_prob = std::min(min_prob, weight[m_yProbIndex]);
            max_prob = std::max(max_prob, weight[m_yProbIndex]);
        }
    }
    // Otherwise we are at a branch point
    else
    {
        std::string    expression = getPredicateExpression(xmlPredicate);
        IXTExprsnNode* xprsnNode  = XprsnParser().parseExpression(expression);
        IXTExprsnNode* ifNode     = parseNode(xmlNodeVec[0]);
        IXTExprsnNode* elseNode   = parseNode(xmlNodeVec[1]);

        //node = new XTIfElseNode<double>(sNodeId, *xprsnNode, *ifNode, *elseNode);
        node = new XTIfElseNode<CTOutPut>(sNodeId, *xprsnNode, *ifNode, *elseNode);
    }

    --nesting_level;

    return node;
}

std::string xmlNewPredictEngineFactory::getPredicateExpression(DOMElement* xmlPredicate)
{
    std::string expression = "";

    if (xmlBase::Dom::getTagName(xmlPredicate) == "CompoundPredicate")
    {
        std::string compoundOperator = xmlBase::Dom::getAttribute(xmlPredicate, "booleanOperator");

        if      (compoundOperator == "or" ) compoundOperator = "|";
        else if (compoundOperator == "and") compoundOperator = "&";
        else
        {
            throw Exception("PredictEngineFactory found unrecognized compound expression: " + compoundOperator);
        }

        std::vector<DOMElement *> xmlPredicateVec;
        xmlBase::Dom::getChildrenByTagName(xmlPredicate, "SimplePredicate", xmlPredicateVec);

        //if (xmlPredicateVec.size() > 2)
        //{
            // Can happen legitimately
            //throw Exception("PredictEngineFactory found simple predicate vector > 2 elements");
        //}

        std::vector<DOMElement *>::iterator xmlPredicateVecIter = xmlPredicateVec.begin();

        expression = getPredicateExpression(*xmlPredicateVecIter++) 
                   + compoundOperator  
                   + getPredicateExpression(*xmlPredicateVecIter++);

        // Pick up any more expressions there might be
        while (xmlPredicateVecIter != xmlPredicateVec.end())
        {
            expression += compoundOperator + getPredicateExpression(*xmlPredicateVecIter++);
        }
    }
    else
    {
        std::string varname   = xmlBase::Dom::getAttribute(xmlPredicate, "field");
        std::string valname   = xmlBase::Dom::getAttribute(xmlPredicate, "value");
        std::string sOperator = xmlBase::Dom::getAttribute(xmlPredicate, "operator");

        // ugliness for now, figure out what to do here...
        if      (sOperator == "greaterOrEqual") sOperator = ">=";
        else if (sOperator == "equal"         ) sOperator = "==";
        else if (sOperator == "lessThan"      ) sOperator = "<";
        else
        {
            throw Exception("PredictEngineFactory found unrecognized operator: " + sOperator);
        }

        //std::string expression = varname + " " + sOperator + " " + valname;
        expression = varname + sOperator + valname;
    }

    return expression;
}


std::string xmlNewPredictEngineFactory::getCTOutputName(const DOMElement* xmlActivityNode)
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
