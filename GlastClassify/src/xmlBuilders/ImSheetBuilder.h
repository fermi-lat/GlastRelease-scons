/** @file Tree.h
*   @brief  Builds DecisionTrees from a given input file.
*
*  $Header$
*/

#ifndef ImSheetBuilder_h
#define ImSheetBuilder_h

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <xercesc/util/XercesDefs.hpp>
#include "../ImActivityNodes/IImActivityNode.h"

XERCES_CPP_NAMESPACE_BEGIN
class  DOMDocument;
class  DOMElement;
XERCES_CPP_NAMESPACE_END

using XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument;
using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;

class DecisionTreeBuilder;
/** @class ImSheetBuilder
*  @brief  Apply a classification Tree derived from Insightful Miner .
*
*  The idea is to encapuslate the xml so that the user only
*   has to provide a filename or string and we will give them access to the
*   linked list of the data which will have the ability to apply the tree 
*   classification algorithms.
*/

class ImSheetBuilder
{
public:
    typedef std::vector<std::string> StringList;

    /** 
        @param searcher call-back to connect with input data at runtime
        @param log  [cout]  ostream for output
        @param iVerbosity [0] -1 no output; 0 errors only; 1 info 2 debug
    */
    ImSheetBuilder(const DOMDocument* document, std::ostream& log=std::cout);
    ~ImSheetBuilder();

    // Return a list of activity nodes
    std::vector<IImActivityNode*> getActivityINodeVec(std::string& type);

    // For output
    void print(std::ostream& out=std::cout) const;

private:

    // Mapping between a node ID and the ImActivityNode object
    typedef std::map<std::string, IImActivityNode*> idToINodeMap;

    // Mapping between the ActivityNode type and a vector of pointers to these nodes
    typedef std::map<std::string, std::vector<IImActivityNode*> > typeToINodeVecMap;

    // Parse the file to find all Activity Nodes
    int findAllActivityNodes(const DOMDocument* document);

    // Parse the link list to connect the Activity Nodes
    int linkActivityNodes(const DOMDocument* document);

    // Vector of Activity Nodes
    std::vector<IImActivityNode*> m_iNodeVec;

    // The "head" ActivityNode
    IImActivityNode*              m_headNode;

    // Map between the node ID and the Activity Node object
    idToINodeMap                  m_idToINodeMap;

    // Map between ActivityNode type and vector of pointers to objects
    typeToINodeVecMap             m_typeToINodeVecMap;

    // Pointer to DecisionTreeBuilder until better idea 
    DecisionTreeBuilder*          m_builder;

    // output related info
    std::ostream&                 m_log;         //! output to this stream
};

#endif // ifdef CLASSIFY_H_
