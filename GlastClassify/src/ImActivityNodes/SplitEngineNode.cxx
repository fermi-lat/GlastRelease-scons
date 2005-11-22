/** @file ImSheetBuilder.cxx
 *    @brief implementation of classification::Tree; declaration and implementation or its private helper classification::ImSheetBuilder::Node
 *
 *    $Header$
 */

#include "SplitEngineNode.h"
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

// Does the "real" work... 
void SplitEngineNode::print(std::ostream& out, int depth) const
{
    // Output our node ID, type and name
    out << indent(depth) << "ID: " << m_id << ", Type: " << m_type << ", Label: " << m_name << std::endl;

    // What is the expression

    // What is the expression
    out << indent(depth) << indent(2) << "Split Expression: ";

    m_xprsnNode->print(out);

    // What do we set depth to?
    depth = m_nodeMap.size() > 1 ? depth + 1 : depth;

    if (m_nodeMap.size() != 2)
    {
        out << indent(depth) << "--> Does not have two nodes! # nodes = " << m_nodeMap.size() << std::endl;
    }

    // Now follow through with all the nodes we point to
    for(IImActivityNodeMap::const_iterator nodeIter = m_nodeMap.begin(); nodeIter != m_nodeMap.end(); nodeIter++)
    {
        nodeIter->second->print(out, depth);
    }

    return;
}

// Does the "real" work... 
void SplitEngineNode::execute()
{
    // If only one node then execute it, otherwise evaluate expression
    if (m_nodeMap.size() == 1)
    {
        m_nodeMap.begin()->second->execute();
    }
    else
    {
        // Evaluate the expression
        bool result = *(reinterpret_cast<const bool*>((*m_xprsnNode)()));

        if (result)
        {
            m_nodeMap[0]->execute();
        }
        else
        {
            m_nodeMap[1]->execute();
        }
    }
    // Now follow through with all the daughter nodes we point to
    //for(IImActivityNodeVec::const_iterator nodeIter = m_nodeVec.begin(); nodeIter != m_nodeVec.end(); nodeIter++)
    //{
    //    (*nodeIter)->execute();
    //}

    return;
}
