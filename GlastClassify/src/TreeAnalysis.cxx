/**@file TreeAnalysis.cxx

@brief implementation of class TreeAnalysis

$Header$
*/

#include "TreeAnalysis.h"
#include "ImActivityNodes/IImActivityNode.h"
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
    m_xtConstants.clear();

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
    try {
    m_headNode->execute();
    } catch(...)
    {
        int j = 0;
    }

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

        // Only dealing with "CTB" variables here
        if (varName.substr(0,3) == "CTB")
        {
            // De-reference to continuous tuple type
            XTcolumnValBase* basePtr = dataIter->second;

            // Null depending on variable type
            if (basePtr->getType() == "continuous")
            {
                XTcolumnVal<REALNUM>* valPtr = dynamic_cast<XTcolumnVal<REALNUM>*>(basePtr);

                valPtr->setDataValue(0.);
                valPtr->clearValidFlag();
            }
            else if (basePtr->getType() == "categorical")
            {
                XTcolumnVal<std::string>* valPtr = dynamic_cast<XTcolumnVal<std::string>*>(basePtr);

                std::string null = "null";
                valPtr->setDataValue(null);
                valPtr->clearValidFlag();
            }
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

                const_cast<GlastClassify::Item*>(m_nTupleMap[dataIter->first])->setDataValue(&result);

                REALNUM test = *(m_nTupleMap[dataIter->first]);

                if (test != result)
                {
                    // what do we do here?
                    int j = 0;
                }
            }
            else if (dataIter->second->getType() == "categorical")
            {
                XTcolumnVal<std::string>* colVal = dynamic_cast<XTcolumnVal<std::string>*>(dataIter->second);
                std::string               result = "";

                if (colVal->dataIsValid()) result = *(*colVal)();

                char* charData = const_cast<char*>(result.c_str());

                const_cast<GlastClassify::Item*>(m_nTupleMap[dataIter->first])->setDataValue(charData);
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

        // Try to look up the variable in the input ntuple. If it exists then we are good to go.
        try
        {
            // Try to look up the variable
            const GlastClassify::Item* item = m_lookup.getItem(varName);

            // If it exists the add to the map
            if (item != 0) m_nTupleMap[varName] = item;
        }
        // If variable doesn't exist then an exception is thrown, we catch it here
        catch (std::invalid_argument&)
        {
            // if this is a CTB variable then we want to add it to the ntuple
            if (varName.substr(0,3) == "CTB")
            {
                if (dataIter->second->getType() == "continuous")
                {
                    // Create a new data object and add to the vector
                    float* newVar = new float;

                    // Initialize it to zero just to make sure
                    *newVar = 0.;

                    // Add it to the ntuple
                    m_lookup.addItem(varName, *newVar);

                    // Now add it to the list
                    m_nTupleMap[varName] = m_lookup.getItem(varName);
                }
                else if (dataIter->second->getType() == "categorical")
                {
                    char* newVar = new char[80];

                    memset(newVar, ' ', 80);

                    m_lookup.addItem(varName, *newVar);

                    m_nTupleMap[varName] = m_lookup.getItem(varName);
                }
            }
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
 //   if (m_headNode) m_headNode->print(out, depth);

    return;
}
