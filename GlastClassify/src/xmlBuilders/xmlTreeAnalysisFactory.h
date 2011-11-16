/**@file xmlTreeAnalysisFactory.h

@brief declaration of class xmlTreeAnalysisFactory
@author T. Burnett

$Header$
*/

#ifndef xmlTreeAnalysisFactory_h
#define xmlTreeAnalysisFactory_h

#include "GlastSvc/GlastClassify/ITupleInterface.h"
#include "../ImActivityNodes/IImActivityNode.h"

#include <vector>
#include <map>
#include <iostream>

#include <xercesc/util/XercesDefs.hpp>
XERCES_CPP_NAMESPACE_BEGIN
class  DOMDocument;
XERCES_CPP_NAMESPACE_END
using XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument;

namespace GlastClassify {

class TreeAnalysis;

/** @class xmlTreeAnalysisFactory
@brief A factory for accessing decision trees

*/
class xmlTreeAnalysisFactory  
{
public:

    /** @brief ctor sets up for tree production 
    @param path file path to a folder containing tree data
    @param lookup Instance of a class supplied to by user, which is called back 
    to find address of each variable
    */
    xmlTreeAnalysisFactory(const std::string& path, ITupleInterface& tuple);

    ~xmlTreeAnalysisFactory() {}

    /** @brief This constructs the desired TreeAnalysis object and returns a pointer
               to it. 
    */
    GlastClassify::TreeAnalysis* buildTreeAnalysis();

    int nodeCount()const{return m_iNodeVec.size();}

private:

    // Mapping between a node ID and the ImActivityNode object
    typedef std::map<std::string, IImActivityNode*> idToINodeMap;

    // Mapping between the ActivityNode type and a vector of pointers to these nodes
    typedef std::map<std::string, std::vector<IImActivityNode*> > typeToINodeVecMap;

    // Parse the file to find all output variables
    int findAllOutputVars(GlastClassify::TreeAnalysis* tree);

    // Parse the file to find all Activity Nodes
    int findAllActivityNodes(GlastClassify::TreeAnalysis* tree);

    // Parse the link list to connect the Activity Nodes
    int linkActivityNodes(GlastClassify::TreeAnalysis* tree);

    // Helps to parse delimited strings
    std::string getNextWord(std::string &sList, int &iEnd);

    // The xml document containing the CT's
    DOMDocument*                     m_domDocument;

    // This looks up the values in the output ntuple
    ITupleInterface&                 m_lookup;

    // Local copy of vector of Activity Nodes
    std::vector<IImActivityNode*>    m_iNodeVec;

    // Map between the node ID and the Activity Node object
    idToINodeMap                     m_idToINodeMap;

    // Map between ActivityNode type and vector of pointers to objects
    typeToINodeVecMap                m_typeToINodeVecMap;

    // output related info
    std::ostream&                    m_log;         //! output to this stream
};


} // namespace

#endif
