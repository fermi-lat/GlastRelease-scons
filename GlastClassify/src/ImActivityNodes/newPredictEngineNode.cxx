/** @file ImSheetBuilder.cxx
 *    @brief implementation of classification::Tree; declaration and implementation or its private helper classification::ImSheetBuilder::Node
 *
 *    $Header$
 */

#include "newPredictEngineNode.h"
#include "src/XT/XTtupleVars.h"
#include "src/XT/XprsnTree.h"

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include <set> 


namespace 
{
    std::string indent(int depth)
    {
        std::string ret("  ");; 
        for( int i =0; i < depth; ++i) ret += "  ";
        return ret;
    }
} // anonomous namespace

newPredictEngineNode::newPredictEngineNode(const std::string& type, const std::string& name, const std::string& id) : 
                                           m_type(type), 
                                           m_name(name), 
                                           m_id(id)  
{
    m_nodeMap.clear();
    m_inputVar.clear();
    m_outputVar.clear();
    m_trees.clear();
}

// Does the "real" work... 
void newPredictEngineNode::print(std::ostream& out, int depth) const
{
    // Output our node ID, type and name
    out << indent(depth) << "ID: " << m_id << ", Type: " << m_type << ", Label: " << m_name << std::endl;

    // What do we set depth to?
    depth = m_nodeMap.size() > 1 ? depth + 1 : depth;

    // Now follow through with all the nodes we point to
    for(IImActivityNodeMap::const_iterator nodeIter = m_nodeMap.begin(); nodeIter != m_nodeMap.end(); nodeIter++)
    {
        nodeIter->second->print(out, depth);
    }

    return;
}

// Does the "real" work... 
void newPredictEngineNode::execute()
{
    // Loop through the forest evaluating the trees
    double weightSum   = 0.;
    double runningWght = 0.;

    //int treeNum = 0;
    std::string predict;

    for(TreePairVector::iterator treeItr = m_trees.begin(); treeItr != m_trees.end(); treeItr++)
    {
        IXTExprsnNode* tree = treeItr->first;

        //double result = *(reinterpret_cast<const double*>((*tree)()));
        const CTOutPut* ctOutPut = reinterpret_cast<const CTOutPut*>((*tree)());

        double result = ctOutPut->getYProb();

        predict = ctOutPut->getScore();

        // For printing during running, leaving here to remember how to do it
        //const XTIfElseNode<double>& ifElseNode = dynamic_cast<const XTIfElseNode<double>& >(*tree);
        //std::cout << "Tree #" << treeNum++ << " evaluates to: " << result << std::endl;
        //ifElseNode.printExp(std::cout);

        weightSum += treeItr->second;
        runningWght += treeItr->second * result;
    }

    //std::cout << "Tree " << m_name << ", running sum=" << runningWght 
    //          << ", weightSum=" << weightSum;

    if (weightSum > 0.) runningWght /= weightSum;

    //std::cout << " gives: " << runningWght << std::endl;
    
    m_xtColumnVal->setDataValue(runningWght);
    m_predict->setDataValue(predict);
    
    // Now follow through with all the daughter nodes we point to
    for(IImActivityNodeMap::const_iterator nodeIter = m_nodeMap.begin(); nodeIter != m_nodeMap.end(); nodeIter++)
    {
        nodeIter->second->execute();
    }

    return;
}
// Override the ostream << operator. Keeping here as an example....
std::ostream& operator<<(std::ostream& stream, const CTOutPut& node)
{
    stream << node.getScore() << "=" << node.getYProb();
    return stream;
}
