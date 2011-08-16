/**@file TreeAnalysis.h

@brief declaration of class TreeAnalysis
@author T. Burnett

$Header$
*/

#ifndef GlastClassify_TreeAnalysis_h
#define GlastClassify_TreeAnalysis_h

#include "GlastClassify/ITupleInterface.h"

class IImActivityNode;

#include <vector>
#include <utility>
#include <map>

#include "XT/XTtupleVars.h"

namespace GlastClassify {

/** @class TreeAnalysis
@brief A factory for accessing decision trees

*/
class TreeAnalysis  
{
public:

    /** @brief ctor sets up for tree production 
    @param path file path to a folder containing tree data
    @param lookup Instance of a class supplied to by user, which is called back 
    to find address of each variable
    */
    TreeAnalysis(ITupleInterface& tuple);

    ~TreeAnalysis();

    /** @brief This will execute the TreeAnalysis nodes with the current event
    */
    void execute();

    /** @brief This will zero the IM analysis output variables
    */
    void zeroCTvals();

    /** @brief This will zero ALL variables in our tuple
    *          ** Should not be used if analysis running on merit tuple! **
    */
    void zeroAllVals();

    /** @brief This stores IM analysis output variables into output ntuple
    */
    void storeCTvals();

    /** @brief Look up the value of a variable stored in the local tuple
               This should be the "standard" method for retrieving information
    */
    double getTupleVal(const std::string& name);

    /** @brief This gives full access to the XTtupleMap
               This is really intended to be used by the TreeAnalysis "Builder"
               not by casual users!
    */
    XTtupleMap& xtTupleMap() {return m_xtTupleMap;}

    /** @brief This gives full access to the XTtupleMap used to hold constants
               This is really intended to be used by the TreeAnalysis "Builder"
               not by casual users!
    */
    XTtupleMap& xtConstants() {return m_xtConstants;}

    /** @brief Used by the TreeAnalysisBuilder to add a new node to our structure
    */
    void addNewNode(IImActivityNode* node);

    /** @brief Used by the TreeAnalysisBuilder to cross reference the ntuple vars
    */
    void crossRefNtupleVars();

    /** @brief Used by the TreeAnalysisBuilder, once all nodes have been added,
               to set the "head" node (the first to be executed
    */
    void setHeadNode(IImActivityNode* head) {m_headNode = head;}

    // For output
    void print(std::ostream& out=std::cout) const;

private:

    // Mapping between the ActivityNode type and a vector of pointers to these nodes
    typedef std::map<std::string, std::vector<IImActivityNode*> > typeToINodeVecMap;

    // The "head" ActivityNode
    IImActivityNode*                   m_headNode;

    // This looks up the values in the output ntuple
    ITupleInterface&                   m_lookup;

    // Class needed to calcluate local variables used in CT's
    XTtupleMap    m_xtTupleMap;

    // Another to hold constants used in CT's
    XTtupleMap    m_xtConstants;

    // Provide map between "local" variables in analysis and 
    // those existing in the input ntuple
    std::map<std::string,const GlastClassify::Item*> m_nTupleMap;

    // Vector of Activity Nodes
    std::vector<IImActivityNode*>    m_iNodeVec;

    // Map between ActivityNode type and vector of pointers to objects
    typeToINodeVecMap                m_typeToINodeVecMap;
};


} // namespace

#endif
