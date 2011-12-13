/**
 * @class CalRawEnergySelectionTool
 *
 * @brief Implements a Gaudi Tool for selecting the best energy reconstruction method 
 *        from the choices available. This is driven off of a particular Tree/CalCluster
 *        association. 
 *
 * @author The Tracking Software Group
 *
 * $Header$
 */

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TreeClusterRelation.h"
#include "Event/Recon/TkrRecon/TkrTree.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/CalRecon/CalEventEnergy.h"

#include "GlastSvc/GlastClassify/IClassifyTool.h"

#include "ICalEnergySelectionTool.h"

class CalRawEnergySelectionTool : public AlgTool, virtual public ICalEnergySelectionTool
{
public:
    /// Standard Gaudi Tool interface constructor
    CalRawEnergySelectionTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~CalRawEnergySelectionTool() {}

    /// @brief Method to fit a single candidate track. Will retrieve any extra info 
    ///        needed from the TDS, then create and use a new KalFitTrack object to 
    ///        fit the track via a Kalman Filter. Successfully fit tracks are then 
    ///        added to the collection in the TDS.
    StatusCode initialize();

    /// @brief Defines a method use to set the energies for the tracks before the
    ///        track fit is run
    const Event::CalCorToolResult* selectBestEnergy(Event::CalEventEnergy*      calEnergy,
                                                    Event::TreeClusterRelation* treeClusRel);


private:
    /// Internal methods
    void   setTupleValues(Event::TkrTree*    tree, 
                          Event::CalCluster* calCluster);
};

//static ToolFactory<CalRawEnergySelectionTool> s_factory;
//const IToolFactory& CalRawEnergySelectionToolFactory = s_factory;
DECLARE_TOOL_FACTORY(CalRawEnergySelectionTool);

//
// Feeds Combo pattern recognition tracks to Kalman Filter
//

CalRawEnergySelectionTool::CalRawEnergySelectionTool(const std::string& type, const std::string& name, const IInterface* parent) :
                 AlgTool(type, name, parent)
{
    //Declare the additional interface
    declareInterface<ICalEnergySelectionTool>(this);

    return;
}

StatusCode CalRawEnergySelectionTool::initialize()
{
    // Purpose and Method: finds the "Global Event Energy" and constrains the
    //                     first two track energies to sum to it.
    // Inputs:  Calorimeter Energy
    // Outputs: Sets the "constrained" energy for all Candidates
    //          Note: This is the energy that will be used for the 
    //          final track fit. 
    // Dependencies: None
    // Restrictions and Caveats:  None.

    //Always believe in success
    StatusCode sc = StatusCode::SUCCESS;
    
    return sc;
}

const Event::CalCorToolResult* CalRawEnergySelectionTool::selectBestEnergy(Event::CalEventEnergy*      calEnergy,
                                                                           Event::TreeClusterRelation* treeClusRel)
{
    // Purpose and Method: Determine the "best" energy from the available choices
    // Inputs:  Tree/Cluster relation and the list of available reconstructed energies
    // Outputs: The object containing the parameters of the "best" energy 
    // Dependencies: None
    // Restrictions and Caveats:  None.
    //Always believe deep down in the success of this venture we are about to embark on
    StatusCode sc = StatusCode::SUCCESS;

    // Define a null pointer to the object we are going to return
    const Event::CalCorToolResult* calResult = 0;

    // Check that we have some objects in place...
    if (calEnergy)
    {
        // We do one thing here: return the raw energy associated to this cluster
        calResult = calEnergy->findLast("CalRawEnergyTool") ;
    }

    return calResult;
}
