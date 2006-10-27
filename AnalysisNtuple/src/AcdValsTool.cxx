
/** @file AcdValsTool.cxx
@brief Calculates the Adc analysis variables
@author Bill Atwood, Leon Rochester

$Header$
*/

#include "ValBase.h"

#include "GaudiKernel/MsgStream.h"   
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/IToolSvc.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"

#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "Event/Recon/AcdRecon/AcdRecon.h"
#include "Event/Digi/AcdDigi.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
// Point used by AcdDigi
//#include "CLHEP/Geometry/Point3D.h"  //<=== check this
#include "CLHEP/Geometry/Transform3D.h"
// Point used by TKR
//#include "geometry/Point.h"  <=== check this

#include <algorithm>
#include <numeric>
#include <map>

/** @class AcdValsTool
@brief Calculates Acd Values
*/

class AcdValsTool : public ValBase
{
public:

    AcdValsTool( const std::string& type, 
        const std::string& name, 
        const IInterface* parent);

    virtual ~AcdValsTool() { }

    StatusCode initialize();

    StatusCode calculate();

    void tkrHitsCount();

private:

    //Global ACDTuple Items
    float ACD_Total_Energy;
    float ACD_Total_Ribbon_Energy;
    float ACD_Tile_Count; 
    float ACD_Ribbon_Count;
    float ACD_ActiveDist;
    float ACD_ActiveDist3D;
    float ACD_ActiveDist3D_Down;
    float ACD_ActiveDist_Energy;
    float ACD_ActiveDist_Energy_Down;
  
    float ACD_Tkr1ActiveDist;
    float ACD_Tkr1ActiveDist_Down;
    float ACD_Tkr1ActiveDist_Energy;
    float ACD_Tkr1ActiveDist_EnergyDown;

    float ACD_VtxActiveDist;
    float ACD_VtxActiveDist_Down;
    float ACD_VtxActiveDist_Energy;
    float ACD_VtxActiveDist_EnergyDown;

    float ACD_GammaDOCA; 
    float ACD_Corner_DOCA;
    float ACD_Tkr1Ribbon_Dist;
    float ACD_TkrRibbon_Dist;
    float ACD_ActDist3DTop;
    float ACD_ActDist3DR0;
    float ACD_ActDist3DR1;
    float ACD_ActDist3DR2;
    float ACD_ActDist3DTop_Down;
    float ACD_ActDist3DR0_Down;
    float ACD_ActDist3DR1_Down;
    float ACD_ActDist3DR2_Down;
    float ACD_tileTopCount;
    float ACD_tileCount0;
    float ACD_tileCount1;
    float ACD_tileCount2;
    float ACD_tileCount3;
    float ACD_ribbon_ActiveDist;
    float ACD_TkrHitsCountTop;
    float ACD_TkrHitsCountRows[4];

    IGlastDetSvc *m_detSvc;
    double m_vetoThresholdMeV;
    double m_tkrHitsCountCut;


};
namespace {

    // predicate to identify top, (row -1) or  side  (row 0-2)
    class acd_row { 
    public:
        acd_row(int row):m_row(row){}
        bool operator() ( std::pair<idents::AcdId ,double> entry){
            return m_row==-1? entry.first.face() == 0 : entry.first.row()==m_row;
        }
        int m_row;
    };
    // used by accumulate to get maximum for a given row
    class acd_max_energy { 
    public:
        acd_max_energy(int row):m_row(row), m_acd_row(row){}
        double operator()(double energy, std::pair<idents::AcdId ,double> entry) {
            return m_acd_row(entry)? std::max(energy, entry.second) : energy;
        }
        int m_row;
        acd_row m_acd_row;
    };
}
// Static factory for instantiation of algtool objects
static ToolFactory<AcdValsTool> s_factory;
const IToolFactory& AcdValsToolFactory = s_factory;

// Standard Constructor
AcdValsTool::AcdValsTool(const std::string& type, 
                         const std::string& name, 
                         const IInterface* parent)
                         : ValBase( type, name, parent )
{    
    // Declare additional interface
    declareInterface<IValsTool>(this); 

    // in mm
    declareProperty("tkrHitsCountCut", m_tkrHitsCountCut=250.0);
}

StatusCode AcdValsTool::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());

    if( ValBase::initialize().isFailure()) return StatusCode::FAILURE;

    setProperties();
    // get the services

    // use the pointer from ValBase

    if( serviceLocator() ) {

        // find GlastDevSvc service
        if (service("GlastDetSvc", m_detSvc, true).isFailure()){
            log << MSG::INFO << "Couldn't find the GlastDetSvc!" << endreq;
            log << MSG::INFO << "Will be unable to calculate ACD_TkrHitsCount" << endreq;
            m_vetoThresholdMeV = 0.4;
        } else {
            StatusCode sc = m_detSvc->getNumericConstByName("acd.vetoThreshold", &m_vetoThresholdMeV);
            if (sc.isFailure()) {
                log << MSG::INFO << "Unable to retrieve threshold, setting the value to 0.4 MeV" << endreq;
                m_vetoThresholdMeV = 0.4;

            }
        }
    } else {
        return StatusCode::FAILURE;
    }

    // load up the map

/** @page anatup_vars
    @section adcvalstool AdCValsTool Variables
    Notes
    - Default Doca/ActiveDistance is -2000.
    - Active distance is negative if a track is outside a tile, 
    positive if inside.
    - For variables called AcdNoXXX, "No" means "Number."

<table>
<tr><th> Variable <th> Type <th> Description					
<tr><td> AcdTotalEnergy 	
<td>F<td>   Total energy deposited in ACD
<tr><td> AcdTileCount 
<td>F<td>   Number of tiles fired
<tr><td> AcdActiveDist3D 	
<td>F<td>   Largest active distance of any track to the edge of any tile 
<tr><td> AcdActDistTileEnergy 
<td>F<td>   The deposited MC energy in the corresponding hit tile 

<tr><td> AcdActiveDist3D_Down	
<td>F<td>   Largest active distance of any track to the edge of any tile, down going side of tracks
<tr><td> AcdActDistTileEnergy_Down 
<td>F<td>   The deposited MC energy in the corresponding hit tile, down going side of tracks

<tr><td>AcdTkr1ActiveDist
<td>F<td>   Largest active distance from track 1 to the edge of any tile
<tr><td>AcdTkr1ActiveDist_Down
<td>F<td>   Largest active distance from track 1 to the edge of any tile, down going side of tracks
<tr><td>AcdTkr1ActDistTileEnergy
<td>F<td>   The deposited MC energy in the corresponding hit tile
<tr><td>AcdTkr1ActDistTileEnergy_Down
<td>F<td>   The deposited MC energy in the corresponding hit tile, down going side of tracks

<tr><td>AcdVtxActiveDist
<td>F<td>   Largest active distance from vertex extrapolation to the edge of any tile
<tr><td>AcdVtxActiveDist_Down
<td>F<td>   Largest active distance from vertex extrapolation to the edge of any tile, down going side of tracks
<tr><td>AcdVtxActDistTileEnergy
<td>F<td>   The deposited MC energy in the corresponding hit tile
<tr><td>AcdVtxActDistTileEnergy_Down
<td>F<td>   The deposited MC energy in the corresponding hit tile, down going side of tracks

<tr><td> AcdGammaDoca 
<td>F<td>   Distance of Gamma to the center of the nearest tile 
<tr><td> AcdCornerDoca 
<td>F<td>   Minimum Distance of Closest Approach of best track to the corner side gaps 
<tr><td> AcdTkrRibbonDist
<td>F<td>   Minimum Distance of Closest Approach of any track to any ribbons that cover gaps
<tr><td> AcdTkr1RibbonDist
<td>F<td>   Minimum Distance of Closest Approach to best track to any ribbons that cover gaps 
<tr><td> AcdActDistTop   
<td>F<td>   Smallest active distance of any track to top tiles 
<tr><td> AcdActDistSideRow[0...3] 	
<td>F<td>   Smallest active distance of any track to tiles in side row [0...3] 
<tr><td> AcdNoTop 
<td>F<td>   hit tile count for top 
<tr><td> AcdNoSideRow[0...3] 
<td>F<td>   Hit Tile count for side row [0...3] 
<tr><td> AcdRibbonActDist   
<td>F<td>   Smallest active distance to any ribbon 
            (considered as a straight line of no thickness) 
<tr><td> AcdTkrHitsCountTop 	
<td>F<td>   Count of the number of TkrClusters within a pre-defined distance 
            (default: 250 mm) of the center of the hit top ACD tiles. 
<tr><td> AcdTkrHitsCountR[0...3] 
<td>F<td>   ditto for ACD tiles in side row [0...3] 
</table>
    */


    addItem("AcdTotalEnergy", &ACD_Total_Energy);
    addItem("AcdRibbonEnergy", &ACD_Total_Ribbon_Energy);
    addItem("AcdRibbonCount", &ACD_Ribbon_Count);
    addItem("AcdTileCount",    &ACD_Tile_Count);

    addItem("AcdActiveDist3D",   &ACD_ActiveDist3D);
    addItem("AcdActDistTileEnergy",   &ACD_ActiveDist_Energy);
    addItem("AcdActiveDist3D_Down", &ACD_ActiveDist3D_Down);
    addItem("AcdActDistTileEnergy_Down", &ACD_ActiveDist_Energy_Down);

    addItem("AcdGammaDoca",    &ACD_GammaDOCA);
    addItem("AcdCornerDoca",    &ACD_Corner_DOCA);
    addItem("AcdTkrRibbonDist",    &ACD_TkrRibbon_Dist);
    addItem("AcdTkr1RibbonDist",    &ACD_Tkr1Ribbon_Dist);

    addItem("AcdTkr1ActiveDist", &ACD_Tkr1ActiveDist);
    addItem("AcdTkr1ActiveDist_Down", &ACD_Tkr1ActiveDist_Down);
    addItem("AcdTkr1ActDistTileEnergy", &ACD_Tkr1ActiveDist_Energy);
    addItem("AcdTkr1ActDistTileEnergy_Down", &ACD_Tkr1ActiveDist_EnergyDown);

    addItem("AcdVtxActiveDist", &ACD_VtxActiveDist);
    addItem("AcdVtxActiveDist_Down", &ACD_VtxActiveDist_Down);
    addItem("AcdVtxActDistTileEnergy", &ACD_VtxActiveDist_Energy);
    addItem("AcdVtxActDistTileEnergy_Down", &ACD_VtxActiveDist_EnergyDown);

    addItem("AcdActDist3DTop",&ACD_ActDist3DTop);
    addItem("AcdActDist3DSideRow0",&ACD_ActDist3DR0);
    addItem("AcdActDist3DSideRow1",&ACD_ActDist3DR1);
    addItem("AcdActDist3DSideRow2",&ACD_ActDist3DR2);

    addItem("AcdActDist3DTop_Down",&ACD_ActDist3DTop_Down);
    addItem("AcdActDist3DSideRow0_Down",&ACD_ActDist3DR0_Down);
    addItem("AcdActDist3DSideRow1_Down",&ACD_ActDist3DR1_Down);
    addItem("AcdActDist3DSideRow2_Down",&ACD_ActDist3DR2_Down);

    addItem("AcdNoTop",        &ACD_tileTopCount);
    addItem("AcdNoSideRow0",   &ACD_tileCount0);
    addItem("AcdNoSideRow1",   &ACD_tileCount1);
    addItem("AcdNoSideRow2",   &ACD_tileCount2);   
    addItem("AcdNoSideRow3",   &ACD_tileCount3);   
    addItem("AcdRibbonActDist", &ACD_ribbon_ActiveDist);

    addItem("AcdTkrHitsCountTop", &ACD_TkrHitsCountTop);
    addItem("AcdTkrHitsCountR0", &ACD_TkrHitsCountRows[0]);
    addItem("AcdTkrHitsCountR1", &ACD_TkrHitsCountRows[1]);
    addItem("AcdTkrHitsCountR2", &ACD_TkrHitsCountRows[2]);
    addItem("AcdTkrHitsCountR3", &ACD_TkrHitsCountRows[3]);

    zeroVals();

    return sc;
}


StatusCode AcdValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;

    tkrHitsCount();

    SmartDataPtr<Event::AcdRecon>           pACD(m_pEventSvc,EventModel::AcdRecon::Event);

    // Recover Track associated info. (not currently used 
    //SmartDataPtr<Event::TkrFitTrackCol>  pTracks(m_pEventSvc,EventModel::TkrRecon::TkrFitTrackCol);
    //SmartDataPtr<Event::TkrVertexCol>     pVerts(m_pEventSvc,EventModel::TkrRecon::TkrVertexCol);
    // Recover pointer to ACD info  

    //Make sure we have valid ACD data
    if (pACD)
    {
        // Make a map relating AcdId to energy in the tile
        std::map<idents::AcdId, double> energyIdMap;

        const std::vector<idents::AcdId> tileIds = pACD->getIdCol();
        const std::vector<double> tileEnergies   = pACD->getEnergyCol();

        int maxTiles = tileIds.size();
        for (int i = 0; i<maxTiles; i++) {
            energyIdMap[tileIds[i]] = tileEnergies[i];
        }

        ACD_Total_Energy  = pACD->getEnergy();
        ACD_Total_Ribbon_Energy = pACD->getRibbonEnergy();
        ACD_Tile_Count    = pACD->getTileCount(); 
        ACD_Ribbon_Count  = pACD->getRibbonCount();

        idents::AcdId tileId = pACD->getMinDocaId();

        tileId            = pACD->getMaxActDist3DId();
        ACD_ActiveDist3D    = pACD->getActiveDist3D();
        ACD_ActiveDist_Energy = energyIdMap[tileId];
        ACD_ActiveDist3D_Down = pACD->getActiveDist3D_Down();
        idents::AcdId tileId_Down = pACD->getMaxActDist3DId_Down();
        ACD_ActiveDist_Energy_Down = energyIdMap[tileId_Down];

        ACD_GammaDOCA     = pACD->getGammaDoca();
        ACD_Corner_DOCA   = pACD->getCornerDoca();
        ACD_ribbon_ActiveDist = pACD->getRibbonActiveDist();

	// loop over AcdGaps & get least distance between track extrapolation & ribbon
	ACD_TkrRibbon_Dist = -2000.;
	ACD_Tkr1Ribbon_Dist = -2000.;
	const Event::AcdTkrGapPocaCol& gaps = pACD->getAcdTkrGapPocaCol();
	for ( Event::AcdTkrGapPocaCol::const_iterator itrGap = gaps.begin(); 
	      itrGap != gaps.end(); itrGap++ ) {
	  // only take the upward going side 
	  if ( (*itrGap)->getArcLength() < 0. ) continue;
	  // only take ribbon gaps	  
	  bool isRibbon(false);
	  switch ( (*itrGap)->getId().gapType() ) {
	  case 1: // X_RibbonSide
	  case 2: // Y_RibbonSide
	  case 3: // Y_RibbonTop
	    isRibbon = true;
	    break;
	  default:
	    break;
	  }
	  if ( !isRibbon ) continue;
	  // only look at first two tracks	  
	  if ( (*itrGap)->trackIndex() == 0 ) {
	    if ( (*itrGap)->getDoca() > ACD_Tkr1Ribbon_Dist ) {
	      ACD_Tkr1Ribbon_Dist = (*itrGap)->getDoca();
	    } 
	  } 
	  if ( (*itrGap)->getDoca() > ACD_TkrRibbon_Dist ) {
	    ACD_TkrRibbon_Dist = (*itrGap)->getDoca();
	  }	  
	}


	ACD_Tkr1ActiveDist = -2000;
	ACD_Tkr1ActiveDist_Energy = 0;
	ACD_Tkr1ActiveDist_Down = -2000;
	ACD_Tkr1ActiveDist_EnergyDown= 0;
	ACD_VtxActiveDist = -2000;
	ACD_VtxActiveDist_Energy= 0;
	ACD_VtxActiveDist_Down = -2000;
	ACD_VtxActiveDist_EnergyDown= 0;
	
	int filledTypeMask(0);


	// loop over AcdTkrHitPoca & get least distance sutff
	const Event::AcdTkrHitPocaCol& pocas = pACD->getAcdTkrHitPocaCol();
	for ( Event::AcdTkrHitPocaCol::const_iterator itrPoca = pocas.begin(); 
	      itrPoca != pocas.end(); itrPoca++ ) {

	  // if already fill all 4 types, break
	  if ( filledTypeMask == 15 ) break;
	  
	  // only take tiles
	  if ( ! (*itrPoca)->getId().tile() ) continue;
	  
	  // only take vertex (-1) and best track (0) 
	  if ( (*itrPoca)->trackIndex() > 0 ) continue;
	  
	  // check up-going v. down going and 
	  int fillType(0);
	  if ( (*itrPoca)->getArcLength() < 0. ) fillType += 1; // down going
	  if ( (*itrPoca)->trackIndex() == -1 ) fillType += 2; // vertices

	  // check to see if we already have an activeDistance of that type
	  int checkFilled = filledTypeMask & ( 1 < fillType );
	  if ( checkFilled != 0 ) continue;
	  checkFilled &= ( 1 < fillType );

	  idents::AcdId theId = (*itrPoca)->getId();

	  // fill the right kind of activeDistance
	  switch ( fillType ) {
	  case 0:
	    ACD_Tkr1ActiveDist = (*itrPoca)->getDoca();
	    ACD_Tkr1ActiveDist_Energy = energyIdMap[theId];
	    break;
	  case 1:
	    ACD_Tkr1ActiveDist_Down  = (*itrPoca)->getDoca();
	    ACD_Tkr1ActiveDist_EnergyDown  = energyIdMap[theId];
	    break;
	  case 2:
	    ACD_VtxActiveDist  = (*itrPoca)->getDoca();
	    ACD_VtxActiveDist_Energy  = energyIdMap[theId];
	    break;
	  case 3:
	    ACD_VtxActiveDist_Down  = (*itrPoca)->getDoca();
	    ACD_VtxActiveDist_EnergyDown = energyIdMap[theId];
	    break;
	  }
	}
    
        const std::vector<double> & adist3D = pACD->getRowActDist3DCol();
	if ( adist3D.size() >= 3 ) {
	  ACD_ActDist3DTop = adist3D[0];
	  ACD_ActDist3DR0 = adist3D[1];
	  ACD_ActDist3DR1 = adist3D[2];
	  ACD_ActDist3DR2 = adist3D[3];
	}

        const std::vector<double> & adist3D_Down = pACD->getRowActDist3DCol_Down();
        if ( adist3D_Down.size() >= 3 ) {
	  ACD_ActDist3DTop_Down = adist3D_Down[0];
	  ACD_ActDist3DR0_Down = adist3D_Down[1];
	  ACD_ActDist3DR1_Down = adist3D_Down[2];
	  ACD_ActDist3DR2_Down = adist3D_Down[3];
	}

        //Code from meritAlg.... 
        // get the map of energy vs tile id: have to construct from two parallel vectors
        float m_acd_tileCount[5];
        const std::vector<double> energies = pACD->getEnergyCol();
        const std::vector<idents::AcdId>& ids = pACD->getIdCol();
        std::vector<double>::const_iterator eit = energies.begin();

        std::map<idents::AcdId, double> emap;
        for( std::vector<idents::AcdId>::const_iterator idit = ids.begin(); 
            idit != ids.end() && eit !=energies.end(); ++idit, ++ eit){
                emap[*idit]=*eit;
            }

        // use acd_row predicate to count number of top tiles
        m_acd_tileCount[0] = std::count_if(emap.begin(), emap.end(), acd_row(-1) );

        // use acd_row predicate to count number of tiles per side row
        if(true)for( int row = 0; row<=3; ++row){ 
            m_acd_tileCount[row+1] = std::count_if(emap.begin(), emap.end(), acd_row(row) );
        }
        ACD_tileTopCount = m_acd_tileCount[0];
        ACD_tileCount0 = m_acd_tileCount[1];
        ACD_tileCount1 = m_acd_tileCount[2];
        ACD_tileCount2 = m_acd_tileCount[3];       
        ACD_tileCount3 = m_acd_tileCount[4];       

    } else {
        return StatusCode::FAILURE;
    }

    return sc;
}

void AcdValsTool::tkrHitsCount() {

    // Purpose and Method:  Count the number of TkrClusters within some distance
    //  of hit ACD tiles.
    // Loops over all AcdDigis (skips those that are not "hit") 
    // Then for each TkrCluster, determine the distance between the center of 
    // The ACD tile and the TkrCluster.  If the distance is below m_tkrHitsCountCut
    // Then count this hit, either as top, R0, R1, R2, R3, depending upon what
    // type of ACD tile this is.

    MsgStream log(msgSvc(), name());
    if (!m_detSvc) return;

    SmartDataPtr<Event::AcdDigiCol> acdDigiCol(m_pEventSvc, EventModel::Digi::AcdDigiCol);
    if (!acdDigiCol) return;

    SmartDataPtr<Event::TkrClusterCol> pClusters(m_pEventSvc,EventModel::TkrRecon::TkrClusterCol);
    if (!pClusters) return;


    Event::AcdDigiCol::const_iterator acdDigiIt;
    for (acdDigiIt = acdDigiCol->begin(); acdDigiIt != acdDigiCol->end(); acdDigiIt++) 
    {
        idents::AcdId acdId = (*acdDigiIt)->getId();
        if (!acdId.tile()) continue;  // skip ribbons

        // Use MC Energy cut for now - until we move to veto discrim
        if ((*acdDigiIt)->getEnergy() < m_vetoThresholdMeV) continue; 
        // Try using the veto discriminator
        // If neither veto discrim from either PMT is high - we should skip
        //if ( (!(*acdDigiIt)->getVeto(Event::AcdDigi::A)) && 
        //     (!(*acdDigiIt)->getVeto(Event::AcdDigi::B)) ) continue;

        idents::VolumeIdentifier volId = (*acdDigiIt)->getVolId();
        std::string str;
        std::vector<double> dim;
        StatusCode sc = m_detSvc->getShapeByID(volId, &str, &dim);
        if ( sc.isFailure() ) {
            log << MSG::WARNING << "Failed to retrieve Shape by Id" << endreq;
            continue;
        }
        HepGeom::Transform3D transform;
        sc = m_detSvc->getTransform3DByID(volId, &transform);
        if (sc.isFailure() ) {
            log << MSG::WARNING << "Failed to get transformation" << endreq;
            continue;
        }

        HepPoint3D center(0., 0., 0.);
        HepPoint3D acdCenter = transform * center;

        // Loop over the clusters and calculate the distance from cluster
        // to ACD tile center
        int numClusters = pClusters->size();
        int icluster;
        for (icluster = 0; icluster < numClusters; icluster++) {
            Event::TkrCluster *clusterTds = (*pClusters)[icluster];
            Point clusterPos = clusterTds->position();
            HepPoint3D clusterP(clusterPos.x(), clusterPos.y(), clusterPos.z());
            double dist = clusterP.distance(acdCenter);
            if (dist < m_tkrHitsCountCut) {
                if (acdId.top()) ++ACD_TkrHitsCountTop;
                if (acdId.side()) {
                    unsigned int row = acdId.row();
                    ++ACD_TkrHitsCountRows[row];
                }
            }
        }
    }

    
    return;
}

