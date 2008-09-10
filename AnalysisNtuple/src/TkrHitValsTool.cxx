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
    int Tkr_Cnv_Lyr_Hits;
    int Tkr_numHitsOnTracks;
    int Tkr_numGhosts;
    int Tkr_numToT255s;
    int Tkr_numGhostsOnTracks;
    int Tkr_numToT255sOnTracks;
    int Tkr_numFlaggedTrackHits;
    int Tkr_numWideClusters;
    int Tkr_Max_controller_hits;
    int Tkr_Fst_Cnv_Lyr;
    int Tkr_NCnv_Lyrs_Hit;
    
    int Tkr_HitsPerLyr[_nLayers];

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

/** @page anatup_vars 
    @section tkrhitvalstool TkrHitValsTool Variables

<table>
<tr><th> Variable <th> Type  <th> Description					
<tr><td> TkrNumHits 	
<td>I<td>   Total number of TKR clusters 
<tr><td> TkrNumHitsOnTracks 	
<td>I<td>   Total number of TKR clusters on tracks
<tr><td> TkrNumGhosts 	
<td>I<td>   Total number of TKR ghost clusters 
<tr><td> TkrNumToT255s 	
<td>I<td>   Total number of clusters with ToT==255 
<tr><td> TkrNumGhostsOnTracks 	
<td>I<td>   Total number of ghost clusters on tracks
<tr><td> TkrNumToT255sOnTracks 	
<td>I<td>   Total number of ToT==255 clusters on tracks
<tr><td> TkrNumFlaggedTrackHits 	
<td>I<td>   Total number of track hits flagged because
            the track had some 255 or ghost hits
<tr><td> TkrNumWideClusters	
<td>I<td>   Number of clusters more than 4 strips wide
<tr><td> TkrFirstLayer
<td>I<td>   First layer containing a cluster 
<tr><td> TkrNumLayersHit
<td>I<td>   Total number of hit layers 
<tr><td> TkrHitsInLyrNN, NN=(00,17)   
<td>I<td>   Number of clusters in (bi)layer NN 
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

    addItem("TkrNumHits",             &Tkr_Cnv_Lyr_Hits);       
    addItem("TkrNumHitsOnTracks",     &Tkr_numHitsOnTracks);
    addItem("TkrFirstLayer",          &Tkr_Fst_Cnv_Lyr);        
    addItem("TkrNumLayersHit",        &Tkr_NCnv_Lyrs_Hit);
    addItem("TkrNumGhosts",           &Tkr_numGhosts);
    addItem("TkrNumToT255s",          &Tkr_numToT255s);
    addItem("TkrNumGhostsOnTracks",   &Tkr_numGhostsOnTracks);
    addItem("TkrNumToT255sOnTracks",  &Tkr_numToT255sOnTracks);
    addItem("TkrNumFlaggedTrackHits", &Tkr_numFlaggedTrackHits);
    addItem("TkrNumWideClusters",     &Tkr_numWideClusters);

    int i;
    char buffer[20];
    for(i=0;i<_nLayers;++i) {
        // vars of form TkrHitsInLyrNN, NN = (00,17), for _nLayers==18
        sprintf(buffer, "TkrHitsInLyr%02i",i);
        addItem(buffer, &Tkr_HitsPerLyr[i]);
    }
            
    zeroVals();
    
    return sc;
}


StatusCode TkrHitValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    // Recover Track associated info. 
    SmartDataPtr<Event::TkrClusterCol>   
        pClusters(m_pEventSvc,EventModel::TkrRecon::TkrClusterCol);
    //Make sure we have valid cluster data

    if (!pClusters) return sc;

    int layerIdx;
    for(layerIdx=0;layerIdx<_nLayers;++layerIdx)
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

        Event::TkrClusterColConItr iter = pClusters->begin();
        for(; iter!=pClusters->end();++iter) {
            bool isGhost;
            bool is255;
            Event::TkrCluster* clust = *iter;
            if(isGhost=clust->isSet(Event::TkrCluster::maskGHOST)) Tkr_numGhosts++;
            if(is255=clust->isSet(Event::TkrCluster::mask255))  Tkr_numToT255s++;
            bool onTrack = clust->hitFlagged();
            if(onTrack) {
                Tkr_numHitsOnTracks++;
                if(isGhost) Tkr_numGhostsOnTracks++;
                if(is255) Tkr_numToT255sOnTracks++;
                if(clust->isSet(Event::TkrCluster::maskSAMETRACK)) Tkr_numFlaggedTrackHits++;
                if(clust->size()>4) Tkr_numWideClusters++;
            }
        }
    
    return sc;
}
