/**@file TreeAnalysis.cxx

@brief implementation of class TreeAnalysis

$Header$
*/

#include "TreeAnalysis.h"
#include "src/ImActivityNodes/IImActivityNode.h"
#include <fstream>
#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <cmath>  // for M_PI, among others

using namespace GlastClassify;

TreeAnalysis::TreeAnalysis(ITupleInterface& tuple)
                         : m_lookup(tuple)
{
    // Clear everything just to be sure
    m_iNodeVec.clear();
    m_typeToINodeVecMap.clear();
    m_xtTupleMap.clear();

    return;
}

TreeAnalysis::~TreeAnalysis()
{
    ///@TODO: delete trees
    return;
}

double TreeAnalysis::getTupleVal(const std::string& name) 
{
    double value = 0;
        
    XTcolumnVal<double>::XTtupleMap::iterator dataIter = m_xtTupleMap.find(name);       
    if (dataIter != m_xtTupleMap.end()) value = *(*(dataIter->second))();

    return value;
}


void TreeAnalysis::execute()
{
    // Transfer tuple variables to local tuple, if needed
    for(XTcolumnVal<double>::XTtupleMap::iterator dataIter = m_xtTupleMap.begin(); 
        dataIter != m_xtTupleMap.end(); dataIter++)
    {
        // Recover the variable name
        const std::string& varName = dataIter->first;

        // Look up cross reference of this variable to the input ntuple
        std::map<std::string,const GlastClassify::Item*>::const_iterator nTupleIter = m_nTupleMap.find(varName);

        // If the cross reference exists, set the local value
        if (nTupleIter != m_nTupleMap.end()) dataIter->second->setDataValue(*(nTupleIter->second));
        // Otherwise, no cross reference... set value to zero and set "valid" flag to false 
        else
        {
            dataIter->second->setDataValue(0.);
            dataIter->second->clearValidFlag();
        }
    }

    // Execute the sheet
    m_headNode->execute();

    // Make a pass through to clear out the "invalid" variables
    // If we are not writing rows then this is probably not an issue...
    for(XTcolumnVal<double>::XTtupleMap::iterator dataIter = m_xtTupleMap.begin(); 
        dataIter != m_xtTupleMap.end(); dataIter++)
    {
        if (!(dataIter->second->dataIsValid())) dataIter->second->setDataValue(0.);
    }

    // Done
    return;
}

/// Add new node to our collection
void TreeAnalysis::addNewNode(IImActivityNode* node)
{
    m_iNodeVec.push_back(node);
    m_typeToINodeVecMap[node->getType()].push_back(node);

    // Set up cross reference to the ntuple
    for(XTcolumnVal<double>::XTtupleMap::iterator dataIter = m_xtTupleMap.begin(); 
        dataIter != m_xtTupleMap.end(); dataIter++)
    {
        const std::string& varName = dataIter->first;

        try
        {
            const GlastClassify::Item* item = m_lookup.getItem(varName);

            if (item != 0) m_nTupleMap[varName] = item;
        }
        catch (std::invalid_argument& arg)
        {
            throw arg;
        }
    }

    return;
}

// For output
void TreeAnalysis::print(std::ostream& out) const
{
    // Start with output of details of the sheet
    out << "***********************" << std::endl;
    out << " Number of Activity nodes: " << m_iNodeVec.size() << std::endl;
    out << " Number of Types found:    " << m_typeToINodeVecMap.size() << std::endl;

    for(typeToINodeVecMap::const_iterator typeIter = m_typeToINodeVecMap.begin(); 
        typeIter != m_typeToINodeVecMap.end(); typeIter++)
    {
        out << "   Type: " << (*typeIter).first << ", Number of instances: " << (*typeIter).second.size() << std::endl;
    }
    out << std::endl;

    // Output the tuple variables this sheet uses
    int numVars = m_xtTupleMap.size();
    out << "Local Tuple Map size: " << numVars << std::endl;

    for(XTcolumnVal<double>::XTtupleMap::const_iterator dataIter = m_xtTupleMap.begin(); 
        dataIter != m_xtTupleMap.end(); dataIter++)
    {
        out << dataIter->first << ",  " << dataIter->second->getName() << std::endl;
    }
    
    out << std::endl;

    int depth = 0;

    // Begin with the head node and output the nodes in their order:
    if (m_headNode) m_headNode->print(out, depth);

    return;
}
