/** @file TkrHitValsTool.cxx
@brief Calculates the Tkr hit analysis variables
@author Bill Atwood, Leon Rochester

$Header$
*/

// Include files

#include "ValBase.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/IToolSvc.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"

#include "TkrUtil/ITkrQueryClustersTool.h"

namespace {
    const int _nLayers = 18;
}

/*! @class TkrHitValsTool
@brief calculates TkrHit values

@authors Bill Atwood, Leon Rochester
*/

class TkrHitValsTool : public ValBase
{
public:
    
    TkrHitValsTool( const std::string& type, 
        const std::string& name, 
        const IInterface* parent);
    
    virtual ~TkrHitValsTool() { }
    
    StatusCode initialize();
    
    StatusCode calculate();
    
private:
    
    //TkrClusters Tuple Items
    float Tkr_Cnv_Lyr_Hits;
    float Tkr_Max_controller_hits;
    float Tkr_Fst_Cnv_Lyr;
    float Tkr_NCnv_Lyrs_Hit;
    
    float Tkr_HitsPerLyr[_nLayers];

    ITkrQueryClustersTool* m_clusTool;
};

// Static factory for instantiation of algtool objects
static ToolFactory<TkrHitValsTool> s_factory;
const IToolFactory& TkrHitValsToolFactory = s_factory;

// Standard Constructor
TkrHitValsTool::TkrHitValsTool(const std::string& type, 
                               const std::string& name, 
                               const IInterface* parent)
                               : ValBase( type, name, parent )
{    
    // Declare additional interface
    declareInterface<IValsTool>(this); 
}

/** @page anatup_vars_optional 
    @section tkrhitvalstool TkrHitValsTool Variables

<table>
<tr><th> Variable <th> Type  <th> Description					
<tr><td> TkrNumHits 	
<td>F<td>   Total number of TKR clusters 
<tr><td> TkrFirstLayer
<td>F<td>   First layer containing a cluster 
<tr><td> TkrNumLayersHit
<td>F<td>   Total number of hit layers 
<tr><td> TkrHitsInLyrNN, NN=(00,17)   
<td>F<td>   Number of clusters in (bi)layer NN 
           (numbered from the bottom of the tracker) 
</table>
*/

StatusCode TkrHitValsTool::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());

    if( ValBase::initialize().isFailure()) return StatusCode::FAILURE;

    // get the services
    
    if( serviceLocator() ) {
        
    } else {
        return StatusCode::FAILURE;
    }

    if ((sc = toolSvc()->retrieveTool("TkrQueryClustersTool", m_clusTool)).isFailure())
    {
        throw GaudiException("Service [TkrQueryClustersTool] not found", name(), sc);
    }
    
    // load up the map

    addItem("TkrNumHits",            &Tkr_Cnv_Lyr_Hits);       
    addItem("TkrFirstLayer",         &Tkr_Fst_Cnv_Lyr);        
    addItem("TkrNumLayersHit",       &Tkr_NCnv_Lyrs_Hit);      
    
    addItem("TkrHitsInLyr00",      &Tkr_HitsPerLyr[0]);      
    addItem("TkrHitsInLyr01",      &Tkr_HitsPerLyr[1]);      
    addItem("TkrHitsInLyr02",      &Tkr_HitsPerLyr[2]);    
    addItem("TkrHitsInLyr03",      &Tkr_HitsPerLyr[3]);      
    addItem("TkrHitsInLyr04",      &Tkr_HitsPerLyr[4]);      
    addItem("TkrHitsInLyr05",      &Tkr_HitsPerLyr[5]);      
    addItem("TkrHitsInLyr06",      &Tkr_HitsPerLyr[6]);      
    addItem("TkrHitsInLyr07",      &Tkr_HitsPerLyr[7]);      
    addItem("TkrHitsInLyr08",      &Tkr_HitsPerLyr[8]);      
    addItem("TkrHitsInLyr09",      &Tkr_HitsPerLyr[9]);      
    addItem("TkrHitsInLyr10",     &Tkr_HitsPerLyr[10]);     
    addItem("TkrHitsInLyr11",     &Tkr_HitsPerLyr[11]);     
    addItem("TkrHitsInLyr12",     &Tkr_HitsPerLyr[12]);     
    addItem("TkrHitsInLyr13",     &Tkr_HitsPerLyr[13]);     
    addItem("TkrHitsInLyr14",     &Tkr_HitsPerLyr[14]);     
    addItem("TkrHitsInLyr15",     &Tkr_HitsPerLyr[15]);     
    addItem("TkrHitsInLyr16",     &Tkr_HitsPerLyr[16]);     
    addItem("TkrHitsInLyr17",     &Tkr_HitsPerLyr[17]); 
    
    zeroVals();
    
    return sc;
}


StatusCode TkrHitValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    // Recover Track associated info. 
    SmartDataPtr<Event::TkrClusterCol>   
        pClusters(m_pEventSvc,EventModel::TkrRecon::TkrClusterCol);

    if (!pClusters) return sc;

    //Make sure we have valid cluster data
    if (pClusters)
    {
        int layerIdx = _nLayers;
        while(layerIdx--)
        {
            int hitCount = m_clusTool->getClusters(idents::TkrId::eMeasureX,layerIdx).size()
                         + m_clusTool->getClusters(idents::TkrId::eMeasureY,layerIdx).size();
            
            if (hitCount > 0)
            {
                Tkr_Fst_Cnv_Lyr    = layerIdx;
                Tkr_NCnv_Lyrs_Hit += 1;
            }
            
            Tkr_HitsPerLyr[layerIdx] = hitCount;
        }
        
        Tkr_Cnv_Lyr_Hits = pClusters->size();
    }
    
    return sc;
}
