/** @file AtwoodTrees.h
    @brief  Declare class AtwoodTrees

$Header$
*/
#ifndef GlastClassify_AtwoodTrees_h
#define GlastClassify_AtwoodTrees_h

#include "GlastSvc/GlastClassify/ITupleInterface.h"

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
    * @param imfile -- the Insightful Minor file created by W. B. Atwood. If empty, use the default

    Uses the tuple object to access current tuple items, and to create new ones.
    */
    AtwoodTrees( ITupleInterface& tuple, std::ostream& log=std::cout, std::string imfile  ="", bool printTreeInfo=false);

    /** run the prediction nodes on the current tuple instance
    */
    bool execute();  

    ~AtwoodTrees();

private:

    // These are variables used by the code
    const Item*   m_TkrNumTracks;
    const Item*   m_CalEnergyRaw  ;
    const Item*   m_CalCsIRLn   ;  
    const Item*   m_obfGamStatus;
    const Item*   m_eventId;
    const Item*   m_run;

    // These are variables to be output to the ntuple 
    // (in alphabetical order)
    float         m_bestEnergyProb;
    float         m_CORE;
    float         m_evtLogEnergyRaw;

    const Item*   m_AcdActiveDist3D;
    const Item*   m_AcdRibbonActDist;
    const Item*   m_AcdCornerDoca;
    const Item*   m_Tkr1SSDVeto;

    int           m_executeTreeCnt;
    int           m_goodVals;
    int           m_caughtVals;

    // Cut selections
    float         m_calEnergyCut;
    float         m_csiRadLenCut;
    float         m_numTracksCut;

    // stream for output (if any)
    std::ostream& m_log;

    // Pointer to the classification tree analysis
    TreeAnalysis* m_treeAnalysis;
 };

} // namespace GlastClassify

#endif
