/** @file AtwoodTrees.h
    @brief  Declare class AtwoodTrees

$Header$
*/
#ifndef GlastClassify_AtwoodTrees_h
#define GlastClassify_AtwoodTrees_h

#include "GlastClassify/ITupleInterface.h"
#include "GlastClassify/ITreeFactory.h"

#include <string>
#include <iostream>
#include <vector>

namespace GlastClassify { 
    
    class ITupleInterface; 

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
    AtwoodTrees( 
        ITupleInterface& tuple, 
        std::ostream&    log       =std::cout, 
        std::string      treepath  =""
        );

    /** run the prediction nodes on the current tuple instance
    */
    void execute();  

    ~AtwoodTrees();

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

    // output quantities: corresponding tuple items
    float m_goodCalProb; 
    float m_goodPsfProb; 
    float m_vtxProb; // vertex or track choice
    float m_gammaProb;
    float m_gammaType;
    float m_BestEnergy;

    std::ostream& m_log;

    GlastClassify::ITreeFactory* m_factory;
 };

} // namespace GlastClassify

#endif
