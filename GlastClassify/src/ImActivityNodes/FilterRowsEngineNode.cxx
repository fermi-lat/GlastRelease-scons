/** @file ImSheetBuilder.cxx
 *    @brief implementation of classification::Tree; declaration and implementation or its private helper classification::ImSheetBuilder::Node
 *
 *    $Header$
 */

#include "FilterRowsEngineNode.h"

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
void FilterRowsEngineNode::print(std::ostream& out, int depth) const
{
    // Output our node ID, type and name
    out << indent(depth) << "ID: " << m_id << ", Type: " << m_type << ", Label: " << m_name << std::endl;

    // What is the expression
    //out << indent(depth+2) << "Filter expression = " << m_expression << std::endl;

    // And again...
    out << indent(depth+2) << "Filter expression: ";
    for(StringList::const_iterator expIter = m_parsedExpression.begin(); expIter != m_parsedExpression.end(); expIter++)
    {
        out << *expIter;
    }
    out << std::endl;

    // What do we set depth to?
    depth = m_nodeVec.size() > 1 ? depth + 1 : depth;

    // Now follow through with all the nodes we point to
    for(IImActivityNodeVec::const_iterator nodeIter = m_nodeVec.begin(); nodeIter != m_nodeVec.end(); nodeIter++)
    {
        (*nodeIter)->print(out, depth);
    }

    return;
}
