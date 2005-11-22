/** @file ImSheetBuilder.cxx
 *    @brief implementation of classification::Tree; declaration and implementation or its private helper classification::ImSheetBuilder::Node
 *
 *    $Header$
 */

#include "CreateColumnsEngineNode.h"
#include "src/XT/XprsnTree.h"
#include "src/XT/XTtupleVars.h"

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

// Does the "real" work... 
void CreateColumnsEngineNode::print(std::ostream& out, int depth) const
{
    // Output our node ID, type and name
    out << indent(depth) << "ID: " << m_id << ", Type: " << m_type << ", Label: " << m_name << std::endl;

    for(XprsnNodeMap::const_iterator nodeIter = m_xprsnNodeMap.begin(); nodeIter != m_xprsnNodeMap.end(); nodeIter++)
    {
        out << indent(depth) << indent(2) << nodeIter->first->getName() << " = ";

        nodeIter->second->print(out);
    }

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
void CreateColumnsEngineNode::execute()
{
    // Iterate over the "action" nodes and execute them
    for(XprsnNodeMap::const_iterator nodeIter = m_xprsnNodeMap.begin(); nodeIter != m_xprsnNodeMap.end(); nodeIter++)
    {
        XTcolumnVal<double>* xtTupleVal = nodeIter->first;
        IXTExprsnNode*       xprsnNode  = nodeIter->second;

        double result = *(reinterpret_cast<const double*>((*xprsnNode)()));

        xtTupleVal->setDataValue(result);

        //std::cout << "Value: " << xtTupleVal->getName() << " = ";
        //xprsnNode->print(std::cout, false);
        //std::cout << " = " << result << std::endl;
        //double res2 = result*result;
    }

    // Now follow through with all the daughter nodes we point to
    for(IImActivityNodeMap::const_iterator nodeIter = m_nodeMap.begin(); nodeIter != m_nodeMap.end(); nodeIter++)
    {
        nodeIter->second->execute();
    }

    return;
}
