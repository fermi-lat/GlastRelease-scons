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

#include <sstream>

using namespace GlastClassify;

TreeAnalysis::TreeAnalysis(ITupleInterface& tuple)
                         : m_lookup(tuple)
{
    // Clear everything just to be sure
    m_iNodeVec.clear();
    m_typeToINodeVecMap.clear();
    m_xtTupleMap.clear();
    m_ctbVarMap.clear();

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
        
    XTtupleMap::iterator dataIter = m_xtTupleMap.find(name);       
    if (dataIter != m_xtTupleMap.end())
    {
        XTcolumnValBase* basePtr = dataIter->second;
        
        if (basePtr->getType() == "continuous") value = *(*(dynamic_cast<XTcolumnVal<REALNUM>*>(basePtr)))();
    }

    return value;
}


void TreeAnalysis::execute()
{
    // Initialize the output variables
    //for(std::map<std::string, float*>::iterator iter = m_ctbVarMap.begin(); iter != m_ctbVarMap.end(); iter++)
    //{
    //    float* ctbVar = iter->second;
    //    *ctbVar = 0;
    //}

    // Transfer tuple variables to local tuple, if needed
    for(XTtupleMap::iterator dataIter = m_xtTupleMap.begin(); dataIter != m_xtTupleMap.end(); dataIter++)
    {
        // Recover the variable name
        const std::string& varName = dataIter->first;

        // Look up cross reference of this variable to the input ntuple
        std::map<std::string,const GlastClassify::Item*>::const_iterator nTupleIter = m_nTupleMap.find(varName);

        // De-reference to continuous tuple type
        XTcolumnValBase* basePtr = dataIter->second;

        if (basePtr->getType() != "continuous") continue;

        XTcolumnVal<REALNUM>* valPtr = dynamic_cast<XTcolumnVal<REALNUM>*>(basePtr);

        // If the cross reference exists, set the local value
        //if (nTupleIter != m_nTupleMap.end()) valPtr->setDataValue(*(nTupleIter->second));
        if (nTupleIter != m_nTupleMap.end() && varName.substr(0,3) != "CTB") 
        {
            // Test precision theory here
            double tempRes = *(nTupleIter->second);

            std::stringstream convString;

            convString.setf(std::ios::fixed);
            convString.precision(8);
            convString << tempRes;

            float tempRes2;

            convString >> tempRes2;

            valPtr->setDataValue(tempRes2);
        }
        // Otherwise, no cross reference... set value to zero and set "valid" flag to false 
        else
        {
            valPtr->setDataValue(0.);
            valPtr->clearValidFlag();
        }
    }

    // Execute the sheet
    m_headNode->execute();

    // Done
    return;
}

// Zero value of output CT variables
void TreeAnalysis::zeroCTvals()
{
    // Transfer tuple variables to local tuple, if needed
    for(XTtupleMap::iterator dataIter = m_xtTupleMap.begin(); dataIter != m_xtTupleMap.end(); dataIter++)
    {
        // Recover the variable name
        const std::string& varName = dataIter->first;

        // De-reference to continuous tuple type
        XTcolumnValBase* basePtr = dataIter->second;

        // Only continuous variables for now...
        if (basePtr->getType() != "continuous") continue;

        XTcolumnVal<REALNUM>* valPtr = dynamic_cast<XTcolumnVal<REALNUM>*>(basePtr);

        // If the cross reference exists, set the local value
        //if (nTupleIter != m_nTupleMap.end()) valPtr->setDataValue(*(nTupleIter->second));
        if (varName.substr(0,3) == "CTB") 
        {
            valPtr->setDataValue(0.);
            valPtr->clearValidFlag();
        }
    }

    return;
}

// Copy CT output variables to the output ntuple
void TreeAnalysis::storeCTvals()
{
    // Copy calculated ctb values to the output ntuple
    for(XTtupleMap::iterator dataIter = m_xtTupleMap.begin(); dataIter != m_xtTupleMap.end(); dataIter++)
    {
        // Check if CTB variable, if so then set data value for tuple
        if (dataIter->first.substr(0,3) == "CTB")
        {
            // Currently only dealing with "continuous" variables
            if (dataIter->second->getType() == "continuous")
            {
                XTcolumnVal<REALNUM>* colVal = dynamic_cast<XTcolumnVal<REALNUM>*>(dataIter->second);
                REALNUM               result = 0.;

                if (colVal->dataIsValid()) result = *(*colVal)();

                float* ctbVarAddr = m_ctbVarMap[dataIter->first];
                
                *ctbVarAddr = result;

                int checkit = 0;
            }
        }
    }

    return;
}


/// Add new node to our collection
void TreeAnalysis::addNewNode(IImActivityNode* node)
{
    m_iNodeVec.push_back(node);
    m_typeToINodeVecMap[node->getType()].push_back(node);

    return;
}

void TreeAnalysis::crossRefNtupleVars()
{
    // Set up cross reference to the ntuple
    for(XTtupleMap::iterator dataIter = m_xtTupleMap.begin(); dataIter != m_xtTupleMap.end(); dataIter++)
    {
        const std::string& varName = dataIter->first;

        // Our first task is to determine if this is going to be output to the tuple
        // The convention for this to be true is that the variable name is prefixed
        // by the letters "CTB"
        if (varName.substr(0,3) == "CTB" && dataIter->second->getType() != "categorical")
        {
            // Create a new data object and add to the vector
            float* newVar = new float;

            // Initialize it to zero just to make sure
            *newVar = 0.;

            m_ctbVarMap[varName] = newVar;

            // Add it to the ntuple
            m_lookup.addItem(varName, *newVar);
        }

        // There seem to be two modes of return here, depending on whether we are
        // using the "real" (ie merit) ntuple or the home brewed on
        // Need a try/catch clause to keep from crashing
        try
        {
            // Try to look up the variable
            const GlastClassify::Item* item = m_lookup.getItem(varName);

            // If it exists the add to the map
            if (item != 0) m_nTupleMap[varName] = item;
        }
        // This happens when using "real" ntuple, if variable doesn't exist
        // an exception is thrown
        catch (std::invalid_argument& arg)
        {
            // No need to do anything
            continue;
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
    out << " " << std::endl;

    for(typeToINodeVecMap::const_iterator typeIter = m_typeToINodeVecMap.begin(); 
        typeIter != m_typeToINodeVecMap.end(); typeIter++)
    {
        out << "   Type: " << (*typeIter).first << ", Number of instances: " << (*typeIter).second.size() << std::endl;
    }
    out << std::endl;

    // Output the tuple variables this sheet uses
    int numVars = m_xtTupleMap.size();
    out << "Local Tuple Map size: " << numVars << std::endl;

    for(XTtupleMap::const_iterator dataIter = m_xtTupleMap.begin(); dataIter != m_xtTupleMap.end(); dataIter++)
    {
        out << dataIter->first << ",  " << dataIter->second->getName();

        // Look up cross reference of this variable to the input ntuple
        std::map<std::string,const GlastClassify::Item*>::const_iterator nTupleIter = m_nTupleMap.find(dataIter->first);

        // If the cross reference exists, set the local value
        if (nTupleIter != m_nTupleMap.end()) out << ", is in the output ntuple";

        out << std::endl;
    }
    
    out << std::endl;

    int depth = 0;

    // Begin with the head node and output the nodes in their order:
    if (m_headNode) m_headNode->print(out, depth);

    return;
}
