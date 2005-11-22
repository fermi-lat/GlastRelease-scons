/** @file AtwoodTrees.h
    @brief  Declare class AtwoodTrees

$Header$
*/
#ifndef GlastClassify_AtwoodTrees_h
#define GlastClassify_AtwoodTrees_h

#include "GlastClassify/ITupleInterface.h"

#include <string>
#include <iostream>
#include <vector>

namespace GlastClassify { 
    
    class ITupleInterface; 
    class TreeAnalysis;

/** @class AtwoodTrees
    @brief Manage Atwood-inspired classification trees, creating new tuple variables
    based on values found in the tuple.
*/
class AtwoodTrees 
{
public:
    /** set up the trees:
    * @param tuple -- abstract interface to a tuple that sets up access to tuple items, and creates new ones
    * @param treepath -- file path to the root of the tree definitions

    Uses the tuple object to access current tuple items, and to create new ones.
    */
    AtwoodTrees( ITupleInterface& tuple, std::ostream& log=std::cout, std::string treepath  ="");

    /** run the prediction nodes on the current tuple instance
    */
    bool execute();  

    ~AtwoodTrees();

private:

    const Item*   m_TkrNumTracks;
    const Item*   m_CalEnergyRaw  ;
    const Item*   m_CalCsIRLn   ;  
    const Item*   m_EvtEventId;

    float         m_bestEnergyProb;
    float         m_profileProb;
    float         m_lastLayerProb;
    float         m_trackerProb;
    float         m_paramProb;
    float         m_CTBestEnergy;
    float         m_CTBdeltaEoE;
    float         m_VTX;
    float         m_CORE;
    float         m_GAM;

    int           m_executeTreeCnt;
    int           m_goodVals;
    int           m_caughtVals;

    // stream for output (if any)
    std::ostream& m_log;

    // Pointer to the classification tree analysis
    TreeAnalysis* m_treeAnalysis;
 };

} // namespace GlastClassify

#endif
