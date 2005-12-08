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

    // These are variables used by the code
    const Item*   m_TkrNumTracks;
    const Item*   m_CalEnergyRaw  ;
    const Item*   m_CalCsIRLn   ;  
    const Item*   m_EvtEventId;

    // These are variables to be output to the ntuple 
    // (in alphabetical order)
    float         m_acdLowerTileCount;
    float         m_acdUpperTileCount;
    float         m_bestPsfErr;
    float         m_bestXDir;
    float         m_bestYDir;
    float         m_bestZDir;
    float         m_bestDeltaEoE;
    float         m_bestEnergy;
    float         m_bestEnergyProb;
    float         m_CORE;
    float         m_calFrontBackRatio;
    float         m_calMaxXtalRatio;
    float         m_evtLogEnergyRaw;
    float         m_GAM;
    float         m_goodEnergy;
    float         m_lastLayerProb;
    float         m_paramProb;
    float         m_profileProb;
    float         m_tkrEnergyFrac;
    float         m_tkrLATEdge;
    float         m_trackerProb;
    float         m_VTX;

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
