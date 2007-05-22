
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

#include "AcdUtil/AcdTileFuncs.h"

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

    void setId(const std::vector<Event::AcdTkrIntersection*> vertexUp, 
                        const std::vector<Event::AcdTkrIntersection*> vertexDown,
                        const std::vector<Event::AcdTkrIntersection*> trackUp,
                        const std::vector<Event::AcdTkrIntersection*> trackDown,
                        unsigned int &retId, bool findRibbon=false);


    void findId(const std::vector<Event::AcdTkrIntersection*> vec, idents::AcdId &retId, bool findRibbon=false);

    void reconId(const Event::AcdRecon *pACD);

private:

    //Global ACDTuple Items
    float ACD_Total_Energy;
    float ACD_Total_Ribbon_Energy;
    float ACD_Tile_Count; 
    float ACD_Ribbon_Count;
    float ACD_ActiveDist3D;
    float ACD_ActiveDist3D_Err;
    float ACD_ActiveDist3D_ArcLen;

    float ACD_ActiveDist3D_Down;
    float ACD_ActiveDist_Energy;
    float ACD_ActiveDist_Energy_Down;
  
    float ACD_Tkr1ActiveDist;
    float ACD_Tkr1ActiveDist_Err;
    float ACD_Tkr1ActiveDist_ArcLen;

    float ACD_Tkr1ActiveDist_Down;
    float ACD_Tkr1ActiveDist_Energy;
    float ACD_Tkr1ActiveDist_EnergyDown;

    float ACD_VtxActiveDist;
    float ACD_VtxActiveDist_ArcLen;

    float ACD_VtxActiveDist_Down;
    float ACD_VtxActiveDist_Energy;
    float ACD_VtxActiveDist_EnergyDown;

    float ACD_Corner_DOCA;
    float ACD_Tkr1Ribbon_Dist;
    float ACD_TkrRibbon_Dist;
    float ACD_Tkr1Hole_Dist;
    float ACD_TkrHole_Dist;
  
    float ACD_ribbon_ActiveDist;
    unsigned int ACD_TileIdRecon;
    unsigned int ACD_RibbonIdRecon;

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
<tr><td> AcdRibbonCount 
<td>F<td>   Number of ribbons fired
<tr><td> AcdActiveDist3D 	
<td>F<td>   Largest active distance of any track to the edge of any tile 
<tr><td> AcdActiveDist3DErr	
<td>F<td>   Error on largest active distance of any track to the edge of any tile 
<tr><td> AcdActiveDist3DArcLen	
<td>F<td>   Arclength from head of track at which active distance was calculated 

<tr><td> AcdActDistTileEnergy 
<td>F<td>   The deposited energy in the corresponding hit tile 

<tr><td> AcdActiveDist3D_Down	
<td>F<td>   Largest active distance of any track to the edge of any tile, down going side of tracks
<tr><td> AcdActDistTileEnergy_Down 
<td>F<td>   The deposited energy in the corresponding hit tile, down going side of tracks

<tr><td>AcdTkr1ActiveDist
<td>F<td>   Largest active distance from track 1 to the edge of any tile
<tr><td>AcdTkr1ActiveDistErr
<td>F<td>   Error on largest active distance from track 1 to the edge of any tile
<tr><td>AcdTkr1ActiveDistArcLen
<td>F<td>   Arclength from head of track 1 at which active distance was calculated for that track
<tr><td>AcdTkr1ActiveDist_Down
<td>F<td>   Largest active distance from track 1 to the edge of any tile, down going side of tracks
<tr><td>AcdTkr1ActDistTileEnergy
<td>F<td>   The deposited energy in the corresponding hit tile
<tr><td>AcdTkr1ActDistTileEnergy_Down
<td>F<td>   The deposited MC energy in the corresponding hit tile, down going side of tracks

<tr><td>AcdVtxActiveDist
<td>F<td>   Largest active distance from vertex extrapolation to the edge of any tile
<tr><td>AcdVtxActiveDistArcLen
<td>F<td>   Arclength along vertex reconstruction at which active distance was calculated for vertex.
<tr><td>AcdVtxActiveDist_Down
<td>F<td>   Largest active distance from vertex extrapolation to the edge of any tile, down going side of tracks
<tr><td>AcdVtxActDistTileEnergy
<td>F<td>   The deposited energy in the corresponding hit tile
<tr><td>AcdVtxActDistTileEnergy_Down
<td>F<td>   The deposited energy in the corresponding hit tile, down going side of tracks

<tr><td> AcdCornerDoca 
<td>F<td>   Minimum Distance of Closest Approach of best track to the corner side gaps 
<tr><td> AcdTkrRibbonDist
<td>F<td>   Minimum Distance of Closest Approach of any track to any ribbons that cover gaps
<tr><td> AcdTkr1RibbonDist
<td>F<td>   Minimum Distance of Closest Approach to best track to any ribbons that cover gaps 
<tr><td> AcdTkrHoleDist
<td>F<td>   Minimum Distance of Closest Approach of any track to any of the tile screw holes
<tr><td> AcdTkr1HoleDist
<td>F<td>   Minimum Distance of Closest Approach to best track to any of the tile screw holes
<tr><td> AcdRibbonActDist   
<td>F<td>   Smallest active distance to any ribbon 
            (considered as a straight line of no thickness) 
<tr><td> AcdTileIdRecon
<td>I<td> Tile identifier that was pierced by the reconstructed track.  
          A value of 899 (N/A) is the default and denotes that no ACD tile was 
		  intersected by a reconstructed track.
<tr><td> AcdRibbonIdRecon
<td>I<td> Ribbon identifier that was pierced by the reconstructed track.  
          A value of 899 (N/A) is the default and denotes that no ACD ribbon was 
		  intersected by a reconstructed track.
</table>
    */


    addItem("AcdTotalEnergy", &ACD_Total_Energy);
    addItem("AcdRibbonEnergy", &ACD_Total_Ribbon_Energy);
    addItem("AcdRibbonCount", &ACD_Ribbon_Count);
    addItem("AcdTileCount",    &ACD_Tile_Count);

    addItem("AcdActiveDist3D",   &ACD_ActiveDist3D);
    addItem("AcdActiveDist3DErr",   &ACD_ActiveDist3D_Err);
    addItem("AcdActiveDist3DArcLen",   &ACD_ActiveDist3D_ArcLen); 
    addItem("AcdActDistTileEnergy",   &ACD_ActiveDist_Energy);
    addItem("AcdActiveDist3D_Down", &ACD_ActiveDist3D_Down);
    addItem("AcdActDistTileEnergy_Down", &ACD_ActiveDist_Energy_Down);

    addItem("AcdCornerDoca",    &ACD_Corner_DOCA);
    addItem("AcdTkrRibbonDist",    &ACD_TkrRibbon_Dist);
    addItem("AcdTkr1RibbonDist",    &ACD_Tkr1Ribbon_Dist);
    addItem("AcdTkrHoleDist",    &ACD_TkrHole_Dist);
    addItem("AcdTkr1HoleDist",    &ACD_Tkr1Hole_Dist);

    addItem("AcdTkr1ActiveDist", &ACD_Tkr1ActiveDist);
    addItem("AcdTkr1ActiveDistErr", &ACD_Tkr1ActiveDist_Err);
    addItem("AcdTkr1ActiveDistArcLen", &ACD_Tkr1ActiveDist_ArcLen);

    addItem("AcdTkr1ActiveDist_Down", &ACD_Tkr1ActiveDist_Down);
    addItem("AcdTkr1ActDistTileEnergy", &ACD_Tkr1ActiveDist_Energy);
    addItem("AcdTkr1ActDistTileEnergy_Down", &ACD_Tkr1ActiveDist_EnergyDown);

    addItem("AcdVtxActiveDist", &ACD_VtxActiveDist);
    addItem("AcdVtxActiveDistArcLen", &ACD_VtxActiveDist_ArcLen);
 
    addItem("AcdVtxActiveDist_Down", &ACD_VtxActiveDist_Down);
    addItem("AcdVtxActDistTileEnergy", &ACD_VtxActiveDist_Energy);
    addItem("AcdVtxActDistTileEnergy_Down", &ACD_VtxActiveDist_EnergyDown);

    addItem("AcdRibbonActDist", &ACD_ribbon_ActiveDist);
    addItem("AcdTileIdRecon", &ACD_TileIdRecon);
    addItem("AcdRibbonIdRecon", &ACD_RibbonIdRecon);

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
        reconId(pACD);


        // Make a map relating AcdId to energy in the tile
        std::map<idents::AcdId, double> energyIdMap;

	const Event::AcdHitCol& hitCol = pACD->getAcdHitCol();
	int nHit = hitCol.size();

	static const float MeVMipTile10 = 1.9;
	static const float MeVMipTile12 = 2.28;
	static const float MeVMipRibbon = 0.5;
	
	ACD_Total_Energy = 0.;
	ACD_Total_Ribbon_Energy = 0.;

	for (int iHit(0); iHit < nHit; iHit++ ){
	  const Event::AcdHit* aHit = hitCol[iHit];
	  const idents::AcdId& id = aHit->getAcdId();
	  if ( id.na() ) continue;
	  float mips = aHit->mips();
	  float MeVMip = id.ribbon() ? MeVMipRibbon : 
	    ( id.top() && id.column() == 2 ) ? MeVMipTile12 : MeVMipTile10;
	  float MeV = mips * MeVMip;

	  if ( id.ribbon() ) {
	    ACD_Total_Ribbon_Energy += MeV;
	  } else {
	    ACD_Total_Energy += MeV;
	  }
	  energyIdMap[id] = MeV;
	}
	
        ACD_Tile_Count    = pACD->getTileCount(); 
        ACD_Ribbon_Count  = pACD->getRibbonCount();

        idents::AcdId tileId = pACD->getMinDocaId();

        tileId            = pACD->getMaxActDist3DId();
        ACD_ActiveDist3D    = pACD->getActiveDist3D();
        ACD_ActiveDist_Energy = energyIdMap[tileId];
        ACD_ActiveDist3D_Down = pACD->getActiveDist3D_Down();
        idents::AcdId tileId_Down = pACD->getMaxActDist3DId_Down();
        ACD_ActiveDist_Energy_Down = energyIdMap[tileId_Down];

        ACD_Corner_DOCA   = pACD->getCornerDoca();
        ACD_ribbon_ActiveDist = pACD->getRibbonActiveDist();

	// loop over AcdGaps & get least distance between track extrapolation & ribbon/ hole
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
	ACD_ActiveDist3D_Err = 0.;
	ACD_ActiveDist3D_ArcLen = 0.;
	
	ACD_Tkr1ActiveDist_ArcLen = 0.;
	ACD_Tkr1ActiveDist = -2000;
	ACD_Tkr1ActiveDist_Err = 0.;
	ACD_Tkr1ActiveDist_Energy = 0;
	ACD_Tkr1ActiveDist_Down = -2000;
	ACD_Tkr1ActiveDist_EnergyDown= 0;
	ACD_VtxActiveDist_ArcLen = 0.;
	ACD_VtxActiveDist = -2000;
	ACD_VtxActiveDist_Energy= 0;
	ACD_VtxActiveDist_Down = -2000;
	ACD_VtxActiveDist_EnergyDown= 0;
	
	ACD_TkrHole_Dist = -2000.;
	ACD_Tkr1Hole_Dist = -2000.;	

	int filledTypeMask(0);
	bool isFirst(true);
	float checkActDist(0.);
	double holeDoca(0.), holeDocaError(0.); 
	double planeError(0.);
	int iHole(-1);

	// loop over AcdTkrHitPoca & get least distance sutff
	const Event::AcdTkrHitPocaCol& pocas = pACD->getAcdTkrHitPocaCol();
	for ( Event::AcdTkrHitPocaCol::const_iterator itrPoca = pocas.begin(); 
	      itrPoca != pocas.end(); itrPoca++ ) {

	  const Event::AcdTkrHitPoca* aPoca = (*itrPoca);

	  AcdTileUtil::tileScrewHoleDoca(aPoca->getId(),aPoca->getActiveX(),aPoca->getActiveY(),
					 aPoca->getLocalXXCov(),aPoca->getLocalYYCov(),aPoca->getLocalXYCov(),
					 holeDoca,holeDocaError,iHole);
	  AcdTileUtil::planeErrorProjection(aPoca->getActiveX(),aPoca->getActiveY(),
					    aPoca->getLocalXXCov(),aPoca->getLocalYYCov(),planeError);  

	  if ( isFirst && aPoca->getArcLength() > 0. && aPoca->trackIndex() >= 0 ) {
	    isFirst = false;
	    checkActDist = aPoca->getDoca();
	    if ( fabs(checkActDist-ACD_ActiveDist3D) < 1e-6) {
	      ACD_ActiveDist3D_ArcLen = aPoca->getArcLength();
	      ACD_ActiveDist3D_Err= planeError;
	      ACD_TkrHole_Dist = holeDoca;
	    } else {
	      MsgStream log(msgSvc(), name());
	      log  << "Mismatch between active distance stored in AcdRecon object and first Poca object " 
		   << checkActDist << ' ' << ACD_ActiveDist3D << endreq;
	    }
	  }

	  // if already fill all 4 types, break
	  if ( filledTypeMask == 15 ) break;
	  
	  // only take tiles
	  if ( ! aPoca->getId().tile() ) continue;
	  
	  // only take vertex (-1) and best track (0) 
	  if ( aPoca->trackIndex() > 0 ) continue;
	  
	  // check up-going v. down going and 
	  int fillType(0);
	  if ( aPoca->getArcLength() < 0. )  fillType += 1; // down going	    
	  if ( aPoca->trackIndex() == -1 ) fillType += 2; // vertices
	

	  // check to see if we already have an activeDistance of that type
	  int checkFilled = filledTypeMask & ( 1 << fillType );
	  if ( checkFilled != 0 ) continue;
	  filledTypeMask |= ( 1 << fillType );

	  idents::AcdId theId = aPoca->getId();

	  // fill the right kind of activeDistance
	  switch ( fillType ) {
	  case 0:	    
	    ACD_Tkr1ActiveDist = aPoca->getDoca();	    
	    ACD_Tkr1ActiveDist_Err = planeError;
	    ACD_Tkr1ActiveDist_Energy = energyIdMap[theId];
	    ACD_Tkr1ActiveDist_ArcLen = aPoca->getArcLength();	    
	    ACD_Tkr1Hole_Dist = holeDoca;
	    break;
	  case 1:
	    ACD_Tkr1ActiveDist_Down  = aPoca->getDoca();
	    ACD_Tkr1ActiveDist_EnergyDown  = energyIdMap[theId];
	    break;
	  case 2:
	    ACD_VtxActiveDist  = aPoca->getDoca();
	    ACD_VtxActiveDist_Energy  = energyIdMap[theId];
	    ACD_VtxActiveDist_ArcLen = aPoca->getArcLength();
	    break;
	  case 3:
	    ACD_VtxActiveDist_Down  = aPoca->getDoca();
	    ACD_VtxActiveDist_EnergyDown = energyIdMap[theId];
	    break;
	  }
	}      

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

    /*  RIP 07-05-21  

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

    */
    return;
}


void AcdValsTool::reconId(const Event::AcdRecon *pACD) {
    std::vector<Event::AcdTkrIntersection*> bestTrackUp, bestTrackDown, vertexUp, vertexDown;

    Event::AcdTkrIntersectionCol::const_iterator itrTrackIntersect; 

    idents::AcdId resetId;
    resetId.na(1);
    ACD_TileIdRecon = resetId.id();
    ACD_RibbonIdRecon = resetId.id();

    // loop over the AcdTrackIntersections
    const Event::AcdTkrIntersectionCol& trackIntersectCol = pACD->getAcdTkrIntersectionCol();
    for ( itrTrackIntersect = trackIntersectCol.begin(); 
        itrTrackIntersect != trackIntersectCol.end(); itrTrackIntersect++ ) {

            if ((*itrTrackIntersect)->getTrackIndex() > 0 ) continue;  // could be vertex (-1) or best track (0)
           
            double arcLen = (*itrTrackIntersect)->getArcLengthToIntersection();

            if (((*itrTrackIntersect)->getTrackIndex() == 0)){ // tracks
                if (arcLen > 0) // upward track
                    bestTrackUp.push_back(*itrTrackIntersect);
                else  // downward track
                    bestTrackDown.push_back(*itrTrackIntersect);
            } else if ( ((*itrTrackIntersect)->getTrackIndex() == -1)) { // vertices 
                if (arcLen > 0)
                    vertexUp.push_back(*itrTrackIntersect);
                else 
                    vertexDown.push_back(*itrTrackIntersect);
            }
        }


        setId(vertexUp, vertexDown, bestTrackUp, bestTrackDown, ACD_TileIdRecon, false);
        setId(vertexUp, vertexDown, bestTrackUp, bestTrackDown, ACD_RibbonIdRecon, true);

        bestTrackUp.clear();
        bestTrackDown.clear();
        vertexUp.clear();
        vertexDown.clear();
}


void AcdValsTool::findId(const std::vector<Event::AcdTkrIntersection*> vec, idents::AcdId &retId, bool findRibbon) {

    //Point tilePos, ribbonPos;
    int tileIndex = -1, ribIndex = -1;

    //idents::AcdId tileId;  //, ribId;
    //tileI d.na(1);
    //ribId.na(1);
    retId.na(1);
    if (vec.size() <= 0) {
        //retId = tileId;
        return;
    }

    //if ( (!findRibbon) && (vec[0]->getTileId().tile()) {
    //    tileIndex = 0;
    //    tileId = vec[0]->getTileId();
    //} else {
    //    ribIndex = 0;
    //    ribId = vec[0]->getTileId();
   // }

    // If there is only one TkrIntesection object, then we can return the one id we have
    if (vec.size() == 1) {
        if ( (!findRibbon) && (vec[0]->getTileId().tile()) )
            retId = vec[0]->getTileId();
        else if ( (findRibbon) && (vec[0]->getTileId().ribbon()) ) 
            retId = vec[0]->getTileId();
        return;
    }

    // otherwise, loop over the remaining TkrIntesection objects
    unsigned int ind=0;
    Event::AcdTkrIntersectionCol::const_iterator itrTrackIntersect; 
    for ( itrTrackIntersect = vec.begin(); 
        itrTrackIntersect != vec.end(); itrTrackIntersect++ ) {
            idents::AcdId id = (*itrTrackIntersect)->getTileId();
            if ( (!findRibbon) && (id.tile()) ) {
                if (tileIndex < 0) {  // haven't seen another tile yet
                    tileIndex = ind;
                    retId = id;
                } else { // there was another tile found already
                    // use Z coordinates to choose if one of the found tiles is a "top" tile
                    // chose the greater Z value
                    if ( (retId.top()) || (id.top()) ) {
                        if (vec[tileIndex]->getGlobalPosition().z() < (*itrTrackIntersect)->getGlobalPosition().z()) {
                            retId = (*itrTrackIntersect)->getTileId();
                            tileIndex = ind;
                        }
                    }   
                    // assume side tiles do not overlap in Gleam, so no worries about handling that case right now
                    // right now we just pick up the tile we find first
                }
            } else if (findRibbon) { // ribbon
                if (ribIndex < 0) { // first ribbon intersection found
                    ribIndex = ind;
                    retId = id;
                } 
                // don't worry about overlapping ribbons, just pick up the first ribbon we find
            }

            ind++;
        }


        return;
}

void AcdValsTool::setId(const std::vector<Event::AcdTkrIntersection*> vertexUp, 
                        const std::vector<Event::AcdTkrIntersection*> vertexDown,
                        const std::vector<Event::AcdTkrIntersection*> bestTrackUp,
                        const std::vector<Event::AcdTkrIntersection*> bestTrackDown,
                        unsigned int &retId, bool findRibbon) {
                            
    idents::AcdId vUpId, vDownId, tUpId, tDownId;
    vUpId.na(1);
    vDownId.na(1);
    tUpId.na(1);
    tDownId.na(1);

    findId(vertexUp, vUpId, findRibbon);
    findId(vertexDown, vDownId, findRibbon);
    findId(bestTrackUp, tUpId, findRibbon);
    findId(bestTrackDown, tDownId, findRibbon);

     // Prefer vertices over tracks
     // up versus down
     bool found = false;
     if (!vUpId.na()) {
         retId = vUpId.id();
         found = true;
      } else if (!tUpId.na()) { 
         retId = tUpId.id();
         found = true;
      }

      // no upward intersections found - check downward
      if (!found) {
          if (!vDownId.na()) 
              retId = vDownId.id();
          else if (!tDownId.na()) 
               retId = tDownId.id();
      }
}
