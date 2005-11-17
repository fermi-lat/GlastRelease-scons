/** @file AtwoodLikeTrees.h
    @brief  Declare class AtwoodLikeTrees

$Header$
*/
#ifndef GlastClassify_AtwoodTrees_h
#define GlastClassify_AtwoodTrees_h

#include "GlastClassify/ITupleInterface.h"
#include "GlastClassify/TreeFactory.h"

#include <string>
#include <iostream>
#include <vector>

namespace GlastClassify { 
    
    class ITupleInterface; 

/** @class AtwoodLikeTrees
    @brief Manage Atwood-inspired classification trees, creating new tuple variables
    based on values found in the tuple.
*/
class AtwoodLikeTrees 
{
public:
    /** set up the trees:
    * @param tuple -- abstract interface to a tuple that sets up access to tuple items, and creates new ones
    * @param treepath -- file path to the root of the tree definitions

    Uses the tuple object to access current tuple items, and to create new ones.
    */
    AtwoodLikeTrees( 
        ITupleInterface& tuple, 
        std::ostream&    log       =std::cout, 
        std::string      treepath  =""
        );

    /** run the prediction nodes on the current tuple instance
    */
    void execute();  

    ~AtwoodLikeTrees();

private:

    //! true if the vertex measurment of the gamma direction is better than one-track
    //! must be called after the execute method.
    bool useVertex()const;

    const Item * m_Tkr1FirstLayer;
    const Item * m_CalEnergyRaw  ;
    const Item * m_CalTotRLn   ;  
    const Item * m_VtxAngle    ;  
    const Item * m_EvtEnergyCorr; 
    const Item * m_CalCfpEnergy;  
    const Item * m_CalLllEnergy;  
    const Item * m_CalTklEnergy;  
    const Item * m_TkrEnergyCorr;
    const Item * m_Tkr1Hits     ;
    const Item * m_Tkr2Hits     ;
    const Item * m_Tkr1X0       ;
    const Item * m_Tkr1Y0       ;
    const Item * m_TkrTotalHits ;
    const Item * m_CalXtalMaxEne;

    // output quantities: corresponding tuple items
    // (double for now in order to replace pointer in tuple)
    double m_goodCalProb; 
    double m_goodPsfProb; 
    double m_vtxProb; // vertex or track choice
    double m_gammaProb;
    double m_gammaType;
    double m_BestEnergy;

    // additional values that need to be set for real atwood trees
    float m_BestLogEnergy;
    float m_EvtLogEnergyRaw;
    float m_TkrEnergyFrac;
    float m_TkrTotalHitsNorm;
    float m_TkrTotalHitRatio;
    float m_BestEnergyProb;
    float m_Tkr12DiffHits;
    float m_PSFEneProbPrd;
    float m_CalMaxXtalRatio;
    float m_TkrLATEdge;

    std::ostream& m_log;

    GlastClassify::TreeFactory* m_factory;


    std::vector<GlastClassify::TreeFactory::Tree> m_trees;
 };

} // namespace GlastClassify

#endif
