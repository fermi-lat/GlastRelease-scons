/** @file McValsTool.cxx
@brief Calculates the Mc analysis variables
@author Bill Atwood, Leon Rochester

$Header$
*/
// Include files

#include "ValBase.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/SmartDataLocator.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/IToolSvc.h"
#include "GaudiKernel/IParticlePropertySvc.h"
#include "GaudiKernel/ParticleProperty.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/MCEvent.h"

// TDS class declarations: input data, and McParticle tree
#include "Event/MonteCarlo/McParticle.h"
#include "Event/MonteCarlo/McIntegratingHit.h"
#include "Event/MonteCarlo/McPositionHit.h"

// Reconstructed Tracks.... 
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "Event/Recon/TkrRecon/TkrTree.h"
#include "Event/Recon/TkrRecon/TkrFilterParams.h"

#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/CalRecon/CalClusterMap.h"

#include "Event/Recon/AcdRecon/AcdTkrHitPoca.h"
#include "Event/Recon/AcdRecon/AcdTkrPoint.h"

#include "Event/Recon/AcdRecon/AcdRecon.h"

/*! @class McValsTool
@brief calculates Monte Carlo values

@authors Bill Atwood, Leon Rochester
*/

class McValsTool : public ValBase
{
public:
    
    McValsTool( const std::string& type, 
        const std::string& name, 
        const IInterface* parent);
    
    virtual ~McValsTool() { }
    
    StatusCode initialize();
    
    StatusCode calculate();

  double GetLengthInBox(double *xbound, double *ybound, double *zbound, double *p, double *v);
    
private:

    //Attempt to calculate the energy exiting the tracker
    double getEnergyExitingTkr(Event::McParticle* mcPart);

    //Function to parse the stuff we get from AcdReconAlg
    void getAcdReconVars();
   
    //Pure MC Tuple Items
    float MC_SourceId;
    char  MC_SourceName[80];
    float MC_NumIncident;
    float MC_Id;
    float MC_Charge;
    float MC_Energy;
    float MC_Energy_1;
    float MC_Energy_2;
    float MC_LogEnergy;
    float MC_EFrac;
    float MC_OpenAngle; 
    float MC_RCAngle;
    float MC_AngleMaxE;
    float MC_AngleMinE;
    float MC_AngleAve;
    float MC_TkrExitEne;

    unsigned int   MC_StatusWord;
    
    float MC_x0;
    float MC_y0;
    float MC_z0;
    
    float MC_xdir;
    float MC_ydir;
    float MC_zdir;

    float MC_xdir1;
    float MC_ydir1;
    float MC_zdir1;

    float MC_xdir2;
    float MC_ydir2;
    float MC_zdir2;
  
    // celestial coordinates now set in McCoordsAlg
    

    //MC - Compared to Recon Items

    // TKR
    //float MC_x_err; 
    //float MC_y_err;
    //float MC_z_err;
    
    //float MC_xdir_err; 
    //float MC_ydir_err;
    //float MC_zdir_err;
    
    float MC_dir_err;
    float MC_dir_errN;
    float MC_dir_errN1;
    float MC_TKR1_dir_err;
    float MC_TKR2_dir_err;

    float MC_EvtDeltaEoE;

    // Tree
    float MC_Tree_dir_err;

    // Tree
    float MC_Tree_match_PosX;
    float MC_Tree_match_PosY;
    float MC_Tree_match_PosZ;
    float MC_Tree_match_DirX;
    float MC_Tree_match_DirY;
    float MC_Tree_match_DirZ;
    float MC_Tree_match_dir_err;
    float MC_Tree_match_track1_err;
    float MC_Tree_match_track2_err;
    float MC_Tree_match_id;

    // Filter, if present
    float MC_Filter_dir_err;

    // Best Cal cluster
    float MC_Cal_match_PosX;
    float MC_Cal_match_PosY;
    float MC_Cal_match_PosZ;
    float MC_Cal_match_DirX;
    float MC_Cal_match_DirY;
    float MC_Cal_match_DirZ;
    float MC_Cal_match_energy;
    float MC_Cal_match_rmsTrans;
    float MC_Cal_match_rmsLong;
    float MC_Cal_match_mcDoca;
    float MC_Cal_match_id;

    // ACD
    float MC_AcdXEnter;
    float MC_AcdYEnter;
    float MC_AcdZEnter;

    float MC_AcdActiveDist3D;
    float MC_AcdActDistTileId;
    float MC_AcdActDistTileEnergy;

    double pbound[3][2];
  float MC_length_tkr;
  float MC_length_tkrgap;
  float MC_length_conv_tkr;
  float MC_length_conv_tkrgap;
  float MC_length_cal;
  float MC_length_calgap;

    // to decode the particle charge
    IParticlePropertySvc* m_ppsvc; 

    IValsTool* m_pEvtTool;

};

// Static factory for instantiation of algtool objects
//static ToolFactory<McValsTool> s_factory;
//const IToolFactory& McValsToolFactory = s_factory;
DECLARE_TOOL_FACTORY(McValsTool);

// Standard Constructor
McValsTool::McValsTool(const std::string& type, 
                       const std::string& name, 
                       const IInterface* parent)
                       : ValBase( type, name, parent )
{    
    // Declare additional interface
    declareInterface<IValsTool>(this); 
}

/** @page anatup_vars 
@section mcvalstool McValsTool Variables

<table>
<tr><th> Variable <th> Type <th> Description
<tr><td> McSourceId 
<td>F<td>   Unique integer associated with each MC source type; 
            from McEvent header replaces Mc_src_Id in merit ntuple 
<tr><td> McSourceName 
<td>S<td>   c-string containing the name of the MC source
<tr><td> McNumIncident
<td>F<td>   Number of incident particles, usually 1
            can be zero to N for test beam
<tr><td> McId 
<td>F<td>   StdHepId of primary (-13 = mu+, 22 = gamma, etc.) 
<tr><td> McCharge 
<td>F<td>   Charge of primary 
<tr><td> McEnergy 
<td>F<td>   Kinetic energy of the generated primary particle 
<tr><td> McLogEnergy 
<td>F<td>   log10(McEnergy) 
<tr><td> McEFrac 
<td>F<td>   Fraction of incident energy in highest-energy daughter 
<tr><td> McOpeningAngle 
<td>F<td>   Actual opening angle between the first and second daughters 
            of the promary as generated, (For a primary photon, 
            these will ordinarily be the electron and positron.) 
<tr><td> McTkrExitEne 
<td>F<td>   Attempt to calculate the total energy <strong>leaving</strong> the tracker volume 
<tr><td> McStatusWord
<td>U<td>   Status bits from G4Generator
<tr><td> Mc[X/Y/Z]0 
<td>F<td>   [x/y/z] coordinate of photon conversion or charged particle origin
<tr><td> Mc[X/Y/Z]Dir 
<td>F<td>   [x/y/z] initial direction cosines of primary particle
<tr><td> Mc[X/Y]Err 
<td>F<td>   REMOVED! [x/y] (found) - [x/y] (Mc) (Mc position taken at the z of the 
            found vertex or first hit)
<tr><td> McZErr 
<td>F<td>   REMOVED! z(<strong>actual</strong> vertex or first hit) - McZ0 
<tr><td> Mc[X/Y/Z]DirErr 
<td>F<td>   REMOVED! [x/y/z]dir (found) - [x/y/z]dir (Mc )
<tr><td> McDirErr 
<td>F<td>   Angle between found direction and Mc direction (radians )
<tr><td> McTkr[1/2]DirErr 
<td>F<td>   Angle between direction of [best/second] track and Mc direction (radians) 
<tr><td> McEvtDeltaEoE
<td>F<td>   fractional error of best from McEnergy
<tr><td> McAcd[X/Y/Z]Enter
<td>F<td>   Position where MC particle enters volume surrounded by ACD
<tr><td> McAcdActiveDist3D
<td>F<td>   Largest active distance from MC particle relative to ACD hit tiles
<tr><td> McAcdActDistTileId
<td>F<td>   ID of tile with the largest active distance
<tr><td> McAcdActDistTileEnergy
<td>F<td>   Energy deposited in tile with the largest active distance
</table>
*/


StatusCode McValsTool::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());

    if( ValBase::initialize().isFailure()) return StatusCode::FAILURE;
     
    if( serviceLocator() ) {
        if( service("ParticlePropertySvc", m_ppsvc, true).isFailure() ) {
            log << MSG::ERROR << "Service [ParticlePropertySvc] not found" << endreq;
            return StatusCode::FAILURE;
        }
    } else {
        return StatusCode::FAILURE;
    }

    IToolSvc* pToolSvc = 0; 
    sc = service("ToolSvc", pToolSvc, true);
    if (!sc.isSuccess ()){
        log << MSG::ERROR << "Can't find ToolSvc, no McEvtDeltaEoE" << endreq;
    } else {
        m_pEvtTool = 0;
        sc = pToolSvc->retrieveTool("EvtValsTool", m_pEvtTool);
        if( sc.isFailure() ) {
            log << MSG::INFO << "Unable to find tool: " "EvtValsTool" << endreq;
            log << "Will carry on anyway, McEvtDeltaEoE will not be calculated" << endreq;
        }
    }

    // load up the map

    addItem("McSourceId",     &MC_SourceId,    true);
    addItem("McSourceName",   MC_SourceName,   true);
    addItem("McNumIncident",  &MC_NumIncident, true);
    addItem("McId",           &MC_Id,          true);  
    addItem("McCharge",       &MC_Charge,      true);
    addItem("McEnergy",       &MC_Energy,      true);  
    addItem("McEnergy1",      &MC_Energy_1,    true);  
    addItem("McEnergy2",      &MC_Energy_2,    true);  
    addItem("McLogEnergy",    &MC_LogEnergy,   true);
    addItem("McEFrac",        &MC_EFrac,       true);
    addItem("McOpenAngle",    &MC_OpenAngle,   true);
    addItem("McReCoilAngle",  &MC_RCAngle,   true);
    addItem("McAngleMaxEne",  &MC_AngleMaxE,   true);
    addItem("McAngleMinEne",  &MC_AngleMinE,   true);
    addItem("McAngleAve",     &MC_AngleAve,    true);
    addItem("McTkrExitEne",   &MC_TkrExitEne,  true);

    // added 5/5/09 LSR
    addItem("McStatusWord",   &MC_StatusWord,  true);

    addItem("McEvtDeltaEoE",  &MC_EvtDeltaEoE);   // moved from EvtValsTool 16-May-2012 LSR
    
    addItem("McX0",           &MC_x0,          true);           
    addItem("McY0",           &MC_y0,          true);           
    addItem("McZ0",           &MC_z0,          true);  
    
    addItem("McXDir",         &MC_xdir,        true );         
    addItem("McYDir",         &MC_ydir,        true );         
    addItem("McZDir",         &MC_zdir,        true );         
  
    addItem("McXDir1",         &MC_xdir1,        true );         
    addItem("McYDir1",         &MC_ydir1,        true );         
    addItem("McZDir1",         &MC_zdir1,        true );  

    addItem("McXDir2",         &MC_xdir2,        true );         
    addItem("McYDir2",         &MC_ydir2,        true );         
    addItem("McZDir2",         &MC_zdir2,        true );         
    // removed 5/5/09 LSR
    //addItem("McXErr",         &MC_x_err);        
    //addItem("McYErr",         &MC_y_err);        
    //addItem("McZErr",         &MC_z_err);        
    
    //addItem("McXDirErr",      &MC_xdir_err);     
    //addItem("McYDirErr",      &MC_ydir_err);     
    //addItem("McZDirErr",      &MC_zdir_err);     
    
    addItem("McDirErr",               &MC_dir_err,               true);      
    addItem("McTkr1DirErr",           &MC_TKR1_dir_err,          true); 
    addItem("McTkr2DirErr",           &MC_TKR2_dir_err,          true); 
    addItem("McDirErrN",              &MC_dir_errN,              true); 
    addItem("McDirErrN1",             &MC_dir_errN1,             true); 

    addItem("McTreeDirErr",           &MC_Tree_dir_err,          true);

    addItem("McBestTreePosX",         &MC_Tree_match_PosX,       true);
    addItem("McBestTreePosY",         &MC_Tree_match_PosY,       true);
    addItem("McBestTreePosZ",         &MC_Tree_match_PosZ,       true);
    addItem("McBestTreeDirX",         &MC_Tree_match_DirX,       true);
    addItem("McBestTreeDirY",         &MC_Tree_match_DirY,       true);
    addItem("McBestTreeDirZ",         &MC_Tree_match_DirZ,       true);
    addItem("McBestTreeDirErr",       &MC_Tree_match_dir_err,    true);
    addItem("McBestTreeTrk1Err",      &MC_Tree_match_track1_err, true);
    addItem("McBestTreeTrk2Err",      &MC_Tree_match_track2_err, true);
    addItem("McBestTreeId",           &MC_Tree_match_id,         true);

    addItem("McFilterDirErr",         &MC_Filter_dir_err,        true);

    addItem("McBestCalPosX",          &MC_Cal_match_PosX,        true);
    addItem("McBestCalPosY",          &MC_Cal_match_PosY,        true);
    addItem("McBestCalPosZ",          &MC_Cal_match_PosZ,        true);
    addItem("McBestCalDirX",          &MC_Cal_match_DirX,        true);
    addItem("McBestCalDirY",          &MC_Cal_match_DirY,        true);
    addItem("McBestCalDirZ",          &MC_Cal_match_DirZ,        true);
    addItem("McBestCalEnergy",        &MC_Cal_match_energy,      true);
    addItem("McBestCalRmsTrans",      &MC_Cal_match_rmsTrans,    true);
    addItem("McBestCalRmsLong",       &MC_Cal_match_rmsLong,     true);
    addItem("McBestCalMcDoca",        &MC_Cal_match_mcDoca,      true);
    addItem("McBestCalId",            &MC_Cal_match_id,          true);

    addItem("McAcdXEnter",            &MC_AcdXEnter,             true);
    addItem("McAcdYEnter",            &MC_AcdYEnter,             true);    
    addItem("McAcdZEnter",            &MC_AcdZEnter,             true);

    addItem("McAcdActiveDist3D",      &MC_AcdActiveDist3D,       true);
    addItem("McAcdActDistTileId",     &MC_AcdActDistTileId,      true);
    addItem("McAcdActDistTileEnergy", &MC_AcdActDistTileEnergy,  true);

  addItem("McLengthInTkr",   &MC_length_tkr);
  addItem("McLengthInTkrGap",   &MC_length_tkrgap);
  addItem("McLengthConvInTkr",   &MC_length_conv_tkr);
  addItem("McLengthConvInTkrGap",   &MC_length_conv_tkrgap);
  addItem("McLengthInCal",   &MC_length_cal);
  addItem("McLengthInCalGap",   &MC_length_calgap);
    
    zeroVals();
    
    return sc;
}

StatusCode McValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    // Recover Track associated info. 
    SmartDataPtr<Event::TkrTreeCol>         pTrees(m_pEventSvc, EventModel::TkrRecon::TkrTreeCol);
    
    // Recover Track associated info. 
    SmartDataPtr<Event::TkrTrackCol>   pTracks(m_pEventSvc,EventModel::TkrRecon::TkrTrackCol); 
    SmartDataPtr<Event::TkrVertexCol>  pVerts(m_pEventSvc,EventModel::TkrRecon::TkrVertexCol);
    SmartDataPtr<Event::TkrFilterParamsCol> pFilterCol(m_pEventSvc, EventModel::TkrRecon::TkrFilterParamsCol);
    // Recover MC Pointer
    SmartDataPtr<Event::McParticleCol> pMcParticle(m_pEventSvc, EventModel::MC::McParticleCol);
    // this is avoid creating the object as a side effect!!
    SmartDataLocator<Event::MCEvent>   pMcEvent(m_pEventSvc, EventModel::MC::Event);

    if(pMcEvent) {
        char temp[2] = "_";
        MC_SourceId = pMcEvent->getSourceId();
        strncpy(MC_SourceName, pMcEvent->getSourceName().c_str(),80);
        if (MC_SourceName=="") strncpy(MC_SourceName, temp, 80);
    }
    
    MC_Energy = -1;
    
    if (pMcParticle) {
        
        Event::McParticleCol::const_iterator pMCPrimary = pMcParticle->begin();
        // Skip the first particle... it's for bookkeeping.
        // The second particle is the first real propagating particle.
        // except in the case of beamtest data!!!
        MC_NumIncident = (*pMCPrimary)->daughterList().size();

        // if there are no incident particles
        if (!MC_NumIncident) return sc;

        // ok, go ahead and call the ACD stuff now
        getAcdReconVars();

        // if there is one incident particle, it's okay to use it as before
        // if there are more than one, I don't know what to do, for now
        //     use the mother particle. That way, at least the energy will
        //     be correct.
        if(MC_NumIncident == 1)
        {
            Event::McParticle::StdHepId hepid= (*pMCPrimary)->particleProperty();
            MC_Id = (double)hepid;
            ParticleProperty* ppty = m_ppsvc->findByStdHepID( hepid );
            
            if (ppty) 
            {
                std::string name = ppty->particle(); 
                MC_Charge = ppty->charge();          
            }

            pMCPrimary++;
        }
        
        HepPoint3D Mc_x0;
        CLHEP::HepLorentzVector Mc_p0;
        // launch point for charged particle; conversion point for neutral
        // Let's try changing this to be the first interaction point for all..
       // Mc_x0 = (MC_Charge==0 ? (*pMCPrimary)->finalPosition() : (*pMCPrimary)->initialPosition());
        Mc_x0 = (*pMCPrimary)->finalPosition(); 
        Mc_p0 = (*pMCPrimary)->initialFourMomentum();

        // there's a method v.m(), but it does something tricky if m2<0
        double mass = sqrt(std::max(Mc_p0.m2(),0.0));
        
        Vector Mc_t0 = Vector(Mc_p0.x(),Mc_p0.y(), Mc_p0.z()).unit();
        
        //Pure MC Tuple Items
        MC_Energy = std::max(Mc_p0.t() - mass, 0.0);
        MC_LogEnergy = log10(MC_Energy);
        
        MC_EvtDeltaEoE = -2.;
        float EvtEnergyCorr;
        if(m_pEvtTool) {
            if(m_pEvtTool->getVal("EvtEnergyCorr", EvtEnergyCorr, NOCALC).isSuccess()){
                if (MC_Energy>0) { 
                    MC_EvtDeltaEoE = (EvtEnergyCorr - MC_Energy)/(MC_Energy);
                }
            } 
        }

        MC_x0     = Mc_x0.x();
        MC_y0     = Mc_x0.y();
        MC_z0     = Mc_x0.z();
        
        MC_xdir   = Mc_t0.x();
        MC_ydir   = Mc_t0.y();
        MC_zdir   = Mc_t0.z();

        // convert to (ra, dec)
        // moved to McCoordsAlg to accomodate Interleave
 
        //Attempt to estimate energy exiting the tracker
        MC_TkrExitEne = getEnergyExitingTkr(*pMCPrimary);
        
        MC_StatusWord = (*pMCPrimary)->statusFlags();

        // If no daughters then nothing happened, e.g. gamma traversed LAT without interacting
        if((*pMCPrimary)->daughterList().size() > 0) 
        {
            SmartRefVector<Event::McParticle> daughters   = (*pMCPrimary)->daughterList();
            SmartRef<Event::McParticle>       pp1         = daughters[0]; 
            std::string                       interaction = pp1->getProcess();

            if(interaction == "conv")  // Its a photon conversion; For comptons "compt" or brems "brem"  
            {
                CLHEP::HepLorentzVector     Mc_p1 = pp1->initialFourMomentum();
                SmartRef<Event::McParticle> pp2   = daughters[1];
                CLHEP::HepLorentzVector     Mc_p2 = pp2->initialFourMomentum();

                double e1 = Mc_p1.t();
                double e2 = Mc_p2.t();

                MC_Energy_1 = e1;
                MC_Energy_2 = e2;

               
                
                Vector Mc_t1 = Vector(Mc_p1.x(),Mc_p1.y(), Mc_p1.z()).unit();
                Vector Mc_t2 = Vector(Mc_p2.x(),Mc_p2.y(), Mc_p2.z()).unit();
                MC_xdir1 = Mc_t1.x();
                MC_ydir1 = Mc_t1.y();
                MC_zdir1 = Mc_t1.z();
                MC_xdir2 = Mc_t2.x();
                MC_ydir2 = Mc_t2.y();
                MC_zdir2 = Mc_t2.z();
                Vector Mc_t12= (Mc_t1 + Mc_t2).unit();
                Vector Mc_Recoil = (e1*Mc_t1 + e2*Mc_t2).unit();
                
                double dot_prod = Mc_t1*Mc_t2;
                if(dot_prod > 1.) dot_prod = 1.;
                MC_OpenAngle = acos(dot_prod);

                dot_prod = Mc_t0*Mc_Recoil;
                if(dot_prod > 1.) dot_prod = 1.;
                MC_RCAngle = acos(dot_prod);

                dot_prod = Mc_t0*Mc_t12;
                if(dot_prod > 1.) dot_prod = 1.;
                MC_AngleAve = acos(dot_prod);

                dot_prod = Mc_t0*Mc_t1;
                if(dot_prod > 1.) dot_prod = 1.; 
                float angle1 = acos(dot_prod);
                dot_prod = Mc_t0*Mc_t2;
                if(dot_prod > 1.) dot_prod = 1.;
                float angle2 = acos(dot_prod);

                MC_EFrac = e1/MC_Energy; 
                if(e1 < e2) {
                    MC_EFrac = e2/MC_Energy;
                    MC_AngleMaxE = angle2;
                    MC_AngleMinE = angle1;
                } else {
                    MC_AngleMaxE = angle1;
                    MC_AngleMinE = angle2;
                }
            }  
        }

        // This should really be a test on vertices, not tracks
        if (pTrees)
        {
            // Ok, need to have some trees to do anything
            if (!pTrees->empty())
            {
                // Recover the "best" tree 
                Event::TkrTree* tree = pTrees->front();

                // Get the axis parameters
                const Event::TkrFilterParams* treeAxis = tree->getAxisParams();

                if (treeAxis)
                {
                    const Vector filterDir = -treeAxis->getEventAxis();

                    MC_Tree_dir_err = acos(std::max(-1., std::min(1., filterDir * Mc_t0)));
                }
            }
        }

        // Try the same with the filter (if it exists)
        if (pFilterCol)
        {
            // Ok, need to have some trees to do anything
            if (!pFilterCol->empty())
            {
                // Recover the "best" tree 
                const Event::TkrFilterParams* filter = pFilterCol->front();

                if (filter)
                {
                    const Vector filterDir = -filter->getEventAxis();

                    MC_Filter_dir_err = acos(std::max(-1., std::min(1., filterDir * Mc_t0)));
                }
            }
        }
        
        //Make sure we have valid reconstructed data
        if (pVerts) 
        {
            // Get the first Vertex - First track of first vertex = Best Track
            if(pVerts->empty()) return sc;
            // Build a map between track and vertex for charged vertices only
            std::map<const Event::TkrTrack*, const Event::TkrVertex*> trackToVertexMap;
            
            // Get the first Vertex - First track of first vertex = Best Track
            Event::TkrVertexConPtr pVtxr = pVerts->begin(); 
            Event::TkrVertex*   gamma = *pVtxr++; 
            Point  x0 = gamma->getPosition();
            Vector t0 = gamma->getDirection();

            // Keep track of track to vertex map
            SmartRefVector<Event::TkrTrack>::const_iterator trackItr = gamma->getTrackIterBegin();

            trackToVertexMap[*trackItr] = gamma;

            if (gamma->getStatusBits() & Event::TkrVertex::TWOTKRVTX) 
                        trackToVertexMap[*(++trackItr)] = gamma;

            // removed 5/5/09
            // Reference position errors at the start of recon track(s)
            //double arc_len = (x0.z()-Mc_x0.z())/Mc_t0.z();
            //HepPoint3D x_start = Mc_x0 + arc_len*Mc_t0;
            //MC_x_err  = x0.x()-x_start.x(); 
            //MC_y_err  = x0.y()-x_start.y();
            //// except for z, use the difference between MC conversion point and start of vertex track
            //MC_z_err  = x0.z()-Mc_x0.z();
            //
            //MC_xdir_err = t0.x()-Mc_t0.x(); 
            //MC_ydir_err = t0.y()-Mc_t0.y();
            //MC_zdir_err = t0.z()-Mc_t0.z();

            bool VTX_set = false;
            for(;pVtxr != pVerts->end(); pVtxr++) 
            {
                Event::TkrVertex* vtxN = *pVtxr; 
                if(vtxN->getStatusBits()& Event::TkrVertex::NEUTRALVTX) 
                {
                    Vector tN = vtxN->getDirection();
                    double acostNtMC = acos(tN*Mc_t0);
                    if(!(VTX_set)) 
                    {
                        MC_dir_errN  = acostNtMC;
                        VTX_set = true;
                    }
                    MC_dir_errN1 = acostNtMC; // Assumes last VTX is 1Tkr Neutral Vtx
                }
                // Not neutral
                else
                {
                    trackItr = vtxN->getTrackIterBegin();

                    trackToVertexMap[*trackItr] = vtxN;

                    if (vtxN->getStatusBits() & Event::TkrVertex::TWOTKRVTX) 
                        trackToVertexMap[*(++trackItr)] = vtxN;
                }
            }   

            
            double cost0tMC = t0*Mc_t0;
            
            MC_dir_err  = acos(cost0tMC);
            
            int nParticles = gamma->getNumTracks(); 
            
            SmartRefVector<Event::TkrTrack>::const_iterator pTrack1 = gamma->getTrackIterBegin();  
            const Event::TkrTrack* track_1 = *pTrack1;
            const Event::TkrTrack* track_2 = 0;
            
            Point  x1 = track_1->front()->getPoint(Event::TkrTrackHit::SMOOTHED);
            Vector t1 = track_1->front()->getDirection(Event::TkrTrackHit::SMOOTHED);
            double cost1tMC = t1*Mc_t0;
            
            MC_TKR1_dir_err  = acos(cost1tMC);
            
            // some confusion here... 
            // we need a better way to find the 2nd track given trees
            // was the 2nd track in the 1st vertex before (not so good!)
            // track_2 already points the 2nd best track
            if(nParticles > 1) 
            {
                track_2   = *(pTrack1 + 1);
                Point  x2 = track_2->front()->getPoint(Event::TkrTrackHit::SMOOTHED);
                Vector t2 = track_2->front()->getDirection(Event::TkrTrackHit::SMOOTHED);
                double cost2tMC = t2*Mc_t0;
                
                MC_TKR2_dir_err  = acos(cost2tMC);
            }

            // Now go back to the trees and find the one which is closest in direction to the 
            // incoming MC direction
            SmartDataPtr<Event::TkrTreeCol> treeCol(m_pEventSvc,"/Event/TkrRecon/TkrTreeCol");

            if (treeCol)
            {
                const Event::TkrTree* bestTree  =  0;
                double                bestAngle = -1.;
                int                   bestId    =  1;
                int                   counter   =  0;

                // Go through tree collection and find best match, determined as closest to the first/best
                // track from the tree
                for(Event::TkrTreeCol::const_iterator treeItr  = treeCol->begin();
                                                      treeItr != treeCol->end();
                                                      treeItr++)
                {
                    const Event::TkrTree* tree   = *treeItr;
                    const Event::TkrTrack* track = tree->front();

                    double cosToMc = track->front()->getDirection(Event::TkrTrackHit::SMOOTHED) * Mc_t0;

                    counter++;

                    if (cosToMc > bestAngle)
                    {
                        bestTree  = tree;
                        bestAngle = cosToMc;
                        bestId    = counter;
                    }
                }

                // Did we make a match?
                if (bestTree)
                {
                    Event::TkrTrackVec::const_iterator trackItr = bestTree->begin();

                    const Event::TkrTrack*  track1 = *trackItr++;
                    const Event::TkrTrack*  track2 = 0;
                    const Event::TkrVertex* vertex = trackToVertexMap[track1];

                    if (trackItr != bestTree->end()) track2 = *trackItr;

                    MC_Tree_match_dir_err    = acos(std::max(-1., std::min(1., vertex->getDirection() * Mc_t0)));
                    MC_Tree_match_track1_err = acos(std::max(-1., 
                                          std::min(1., track1->front()->getDirection(Event::TkrTrackHit::SMOOTHED) * Mc_t0)));

                    if (track2)  MC_Tree_match_track2_err = acos(std::max(-1., 
                                          std::min(1., track2->front()->getDirection(Event::TkrTrackHit::SMOOTHED) * Mc_t0)));

                    MC_Tree_match_id = bestId; //bestTree->getHeadNode()->getTreeId();

                    // Get the axis parameters
                    const Event::TkrFilterParams* treeAxis = bestTree->getAxisParams();

                    if (treeAxis)
                    {
                        const Point filterPos  = treeAxis->getEventPosition();
                        const Vector filterDir = treeAxis->getEventAxis();

                        MC_Tree_match_PosX = filterPos.x();
                        MC_Tree_match_PosY = filterPos.y();
                        MC_Tree_match_PosZ = filterPos.z();
                                           
                        MC_Tree_match_DirX = filterDir.x();
                        MC_Tree_match_DirY = filterDir.y();
                        MC_Tree_match_DirZ = filterDir.z();
                    }
                }
                else
                {
                    int probleminrivercity = 0;
                }
            }

            // Now do the same thing trying to find the "best" cluster in the cal cluster collection
            // Recover pointers to CalClusters and Xtals
            SmartDataPtr<Event::CalClusterMap> pCalClusterMap(m_pEventSvc,EventModel::CalRecon::CalClusterMap); 
            Event::CalClusterVec rawClusterVec;
            if(pCalClusterMap) rawClusterVec = (*pCalClusterMap).get(EventModel::CalRecon::CalRawClusterVec);

            if (rawClusterVec.size()>0)
            {
                Event::CalCluster* bestCluster = rawClusterVec.front();
                int                bestId      = 1;
                double             bestDoca    = 100000.;

                // This is a bit tricky since the photons can come from any direction. 
                // Constrain ourselves to consider only photons incident in the upper half hemisphere
                if (Mc_t0.z() < -0.05)
                {
                    // don't forget that the current Mc_t0 is opposite what we need here
                    Mc_t0 = -Mc_t0;

                    // Convert the MC position to a Point
                    Point mcPos(Mc_x0.x(), Mc_x0.y(), Mc_x0.z());

                    // Count our way through the clusters
                    int clusCounter = 0;

                    // Set start and stop iterators 
                    Event::CalClusterVec::iterator calClusIter = rawClusterVec.begin();

                    // Ok, we are going to loop on clusters and search for the cluster which 
                    // has the best distance of closest approach of the MC axis to the centroid
                    while(calClusIter != rawClusterVec.end())
                      {
                        Event::CalCluster* cluster = *calClusIter;

                        // First get the point on the MC axis at the plane of the cluster centroid
                        double arcLenToCalZ = (cluster->getMomParams().getCentroid().z() - Mc_x0.z()) / Mc_t0.z();
                        Point  mcInCalPlane = mcPos + arcLenToCalZ * Mc_t0;

                        // Ok, now get the vector from this point to the centroid
                        Vector mcToCalCent  = cluster->getMomParams().getCentroid() - mcInCalPlane;

                        // Cross product to get distance between the MC trajectory and the cluster
                        Vector docaVec      = mcToCalCent.cross(Mc_t0);

                        // 3-D DOCA
                        double doca         = docaVec.mag();

                        clusCounter++;

                        if (doca < bestDoca)
                        {
                            bestCluster = cluster;
                            bestDoca    = doca;
                            bestId      = clusCounter;
                        }
                        //
                        calClusIter++;
                      }
                }

                if (bestCluster)
                {
                    MC_Cal_match_PosX     = bestCluster->getMomParams().getCentroid().x();
                    MC_Cal_match_PosY     = bestCluster->getMomParams().getCentroid().y();
                    MC_Cal_match_PosZ     = bestCluster->getMomParams().getCentroid().z();
                    MC_Cal_match_DirX     = bestCluster->getMomParams().getAxis().x();
                    MC_Cal_match_DirY     = bestCluster->getMomParams().getAxis().y();
                    MC_Cal_match_DirZ     = bestCluster->getMomParams().getAxis().z();
                    MC_Cal_match_energy   = bestCluster->getMomParams().getEnergy();
                    MC_Cal_match_rmsTrans = bestCluster->getMomParams().getTransRms();
                    MC_Cal_match_rmsLong  = bestCluster->getMomParams().getLongRms();
                    MC_Cal_match_mcDoca   = bestDoca;
                    MC_Cal_match_id       = bestId;

                }
            }
        } 
    }

    int i,j;
    double Xbound[2];
    double Ybound[2];
    double Zbound[2];

    double towerpitch = 374.5;
    double tkrbottom = 25;
    double tkrtop = 625;
    double calbottom = -47.395-8*21.35;
    double caltop = -47.395;
    
    double tkr_gap = 8;
    double cal_gap = 30;
    double p[3];
    p[0] = MC_x0;
    p[1] = MC_y0;
    p[2] = MC_z0;
    double v[3];
    v[0] = MC_xdir;
    v[1] = MC_ydir;
    v[2] = MC_zdir;
    
    double xcenter,ycenter;
    
    Zbound[0] = tkrbottom;
    Zbound[1] = tkrtop;
    double towerhalfwidth = towerpitch/2-tkr_gap;
    MC_length_tkr = 0;
    for(i=0;i<4;++i)
      {
        xcenter = -1.5*towerpitch+towerpitch*i;
        Xbound[0] = xcenter-towerhalfwidth;
        Xbound[1] = xcenter+towerhalfwidth;
        for(j=0;j<4;++j)
          {
            ycenter = -1.5*towerpitch+towerpitch*j;
            Ybound[0] = ycenter-towerhalfwidth;
            Ybound[1] = ycenter+towerhalfwidth;
            MC_length_tkr += GetLengthInBox(Xbound,Ybound,Zbound,p,v);
          }
      }
    Xbound[0] = -2*towerpitch+tkr_gap;
    Xbound[1] = 2*towerpitch-tkr_gap;
    Ybound[0] = -2*towerpitch+tkr_gap;
    Ybound[1] = 2*towerpitch-tkr_gap;
    MC_length_tkrgap = GetLengthInBox(Xbound,Ybound,Zbound,p,v)-MC_length_tkr;
    //
    MC_length_conv_tkr = 0;
    MC_length_conv_tkrgap = 0;
    if(MC_z0>tkrbottom)
      {
        Zbound[0] = tkrbottom;
        Zbound[1] = MC_z0;
        towerhalfwidth = towerpitch/2-tkr_gap;
        MC_length_conv_tkr = 0;
        for(i=0;i<4;++i)
          {
            xcenter = -1.5*towerpitch+towerpitch*i;
            Xbound[0] = xcenter-towerhalfwidth;
            Xbound[1] = xcenter+towerhalfwidth;
            for(j=0;j<4;++j)
              {
                ycenter = -1.5*towerpitch+towerpitch*j;
                Ybound[0] = ycenter-towerhalfwidth;
                Ybound[1] = ycenter+towerhalfwidth;
                MC_length_conv_tkr += GetLengthInBox(Xbound,Ybound,Zbound,p,v);
              }
          }
        Xbound[0] = -2*towerpitch+tkr_gap;
        Xbound[1] = 2*towerpitch-tkr_gap;
        Ybound[0] = -2*towerpitch+tkr_gap;
        Ybound[1] = 2*towerpitch-tkr_gap;
        MC_length_conv_tkrgap = GetLengthInBox(Xbound,Ybound,Zbound,p,v)-MC_length_conv_tkr;
      }
    //
    Zbound[0] = calbottom;
    Zbound[1] = caltop;
    towerhalfwidth = towerpitch/2-cal_gap;
    MC_length_cal = 0;
    for(i=0;i<4;++i)
      {
        xcenter = -1.5*towerpitch+towerpitch*i;
        Xbound[0] = xcenter-towerhalfwidth;
        Xbound[1] = xcenter+towerhalfwidth;
        for(j=0;j<4;++j)
          {
            ycenter = -1.5*towerpitch+towerpitch*j;
            Ybound[0] = ycenter-towerhalfwidth;
            Ybound[1] = ycenter+towerhalfwidth;
            MC_length_cal += GetLengthInBox(Xbound,Ybound,Zbound,p,v);
          }
      }
    //
    Xbound[0] = -2*towerpitch+cal_gap;
    Xbound[1] = 2*towerpitch-cal_gap;
    Ybound[0] = -2*towerpitch+cal_gap;
    Ybound[1] = 2*towerpitch-cal_gap;
    MC_length_calgap = GetLengthInBox(Xbound,Ybound,Zbound,p,v)-MC_length_cal;
    
    return sc;
}

double McValsTool::getEnergyExitingTkr(Event::McParticle* mcPart)
{
    //This function attempts to estimate the energy which exist the tracker during 
    //an event. The idea is to recursively traverse the McParticle tree and look for
    //particles which start in the tracker (so, for example, the electron and positron
    //from a gamma conversion) and then exit, keep track of the energy of these particles.
    //Due to the vagaries of our simulation, several interesting problems can arise and
    //an attempt is made to deal with these (see notes below).
    //A flaw here is that NO attempt is made to account for continuous energy loss of 
    //charged particles in the tracker. Hopefully, a more sophisticated routine will 
    //appear soon.
    double tkrExitEne = 0.;
    double partEnergy = 0.;
    double partLostE  = 0.;

    //idents::VolumeIdentifier initial = mcPart->getInitialId();
    idents::VolumeIdentifier final   = mcPart->getFinalId();

    //This should check that the initial point of the particle is in the tracker
    //and that its final point is outside. So, it started in the tracker and 
    //carried its energy outside (could be anywhere)
    ////if (initial.size() == 9 && initial[0] == 0 && initial[3] == 1 &&
    ////    final.size() == 9 && (final[0] != 0 || final[3] != 1)       )
    if (final.size() == 9 && (final[0] != 0 || final[3] != 1) )
    {
        //Set total initial energy of this particle...
        partEnergy = mcPart->initialFourMomentum().t();
    }

    //Next step is to go through the daughter list to find the tracker exiting
    //energy due to them
    SmartRefVector<Event::McParticle> daughters = mcPart->daughterList();

    for(SmartRefVector<Event::McParticle>::iterator partIter = daughters.begin(); 
                                                    partIter < daughters.end(); 
                                                    partIter++)
    {
        Event::McParticle* daughter = *partIter;

        //First we check to see if the daughter was created in the tracker. 
        idents::VolumeIdentifier daughterInitial = daughter->getInitialId();

        if (daughterInitial.size() == 9 && daughterInitial[0] == 0 && daughterInitial[3] ==1)
        {
            //Keep track of the energy this daughter takes from its parent
            partLostE += daughter->initialFourMomentum().t();

            //Now get the energy exiting the tracker due to this particle
            tkrExitEne += getEnergyExitingTkr(daughter);
        }
    }

    //Ok, now correct the mcPart energy for that it lost traversing the tracker by 
    //creating daughter particles
    partEnergy -= partLostE;

    if (partEnergy < 0.) partEnergy = 0.;

    //Finally, add this to the energy exiting the tracker
    tkrExitEne += partEnergy;

    //This should now be the energy exiting the tracker due to those particles created 
    //in the tracker. 
    return tkrExitEne;
}


void McValsTool::getAcdReconVars() {

  MsgStream log(msgSvc(), name());  

  double bestActDist(-2000.);
  idents::AcdId bestId;
  std::map<idents::AcdId, double> energyIdMap;

  SmartDataPtr<Event::AcdRecon>           pACD(m_pEventSvc,EventModel::AcdRecon::Event);
  if (pACD) {
    // Make a map relating AcdId to energy in the tile
    const std::vector<idents::AcdId>& tileIds = pACD->getIdCol();
    const std::vector<double>& tileEnergies   = pACD->getEnergyCol();
    
    int maxTiles = tileIds.size();
    for (int i = 0; i<maxTiles; i++) {
      energyIdMap[tileIds[i]] = tileEnergies[i];
    }
  }

  // Here we will loop over calculated poca's to find best (largest) active distance
  SmartDataPtr<Event::AcdTkrHitPocaCol> acdTkrHits(m_pEventSvc,EventModel::MC::McAcdTkrHitPocaCol); 
  if ( ! acdTkrHits ) {
    // no poca's found, set guard values.
    log << "no AcdTkrHitPocas found on TDS" << endreq;
    MC_AcdActiveDist3D = bestActDist;
    MC_AcdActDistTileId = bestId.id();
    MC_AcdActDistTileEnergy = -1.;
  } else {  
    for ( Event::AcdTkrHitPocaCol::const_iterator itr = acdTkrHits->begin();
      itr != acdTkrHits->end(); itr++ ) {
      const Event::AcdTkrHitPoca* aHitPoca = *itr;
      if ( aHitPoca->getDoca() > bestActDist ) {
    bestId = aHitPoca->getId();
    bestActDist = aHitPoca->getDoca();
      }
    }  
    // latch values
    MC_AcdActiveDist3D = bestActDist;
    MC_AcdActDistTileId = bestId.id();
    MC_AcdActDistTileEnergy = energyIdMap[bestId];
  }

  // Here we will get the point the MC parent enters the ACD volume
  SmartDataPtr<Event::AcdTkrPointCol> acdPoints(m_pEventSvc,EventModel::MC::McAcdTkrPointCol);
  if ( ! acdPoints || acdPoints->size() < 1 ) {
    // no point's found, set guard values.
    log << "no AcdTkrPoints found on TDS" << endreq;
    MC_AcdXEnter = MC_AcdYEnter = MC_AcdZEnter = -2000.;
  } else {
    const Event::AcdTkrPoint* entryPoint = (*acdPoints)[0];
    const HepPoint3D& thePoint = entryPoint->getGlobalPosition();
    
    // latch values
    MC_AcdXEnter = thePoint.x(); 
    MC_AcdYEnter = thePoint.y(); 
    MC_AcdZEnter = thePoint.z(); 
  }

}

double McValsTool::GetLengthInBox(double *xbound, double *ybound, double *zbound, double *p, double *v)
{
  int i,j,k;
  double pp[3];
  double ppp[2][3];
  double lambda;

  pbound[0][0] = xbound[0];
  pbound[0][1] = xbound[1];
  pbound[1][0] = ybound[0];
  pbound[1][1] = ybound[1];
  pbound[2][0] = zbound[0];
  pbound[2][1] = zbound[1];

  int inside;
  int ii = 0;

  for(i=0;i<3;++i)
    {
      if(v[i]==0) continue;
      //
      for(j=0;j<2;++j)
        {
          lambda = (pbound[i][j]-p[i])/v[i];
          inside = 0;
          for(k=0;k<3;++k)
            {
              pp[k] = p[k]+lambda*v[k];
              if(pp[k]>=pbound[k][0]-1e-6&&pp[k]<=pbound[k][1]+1e-6) ++inside;
            }
          if(inside<3) continue;
          //
          if(ii>0 && fabs(pp[0]-ppp[0][0])<1e-6 && fabs(pp[1]-ppp[0][1])<1e-6 && fabs(pp[2]-ppp[0][2])<1e-6) continue;
          for(k=0;k<3;++k) ppp[ii][k] = pp[k];
          ++ii;
          if(ii==2) break;
        }
      if(ii==2) break;
    }
  if(ii<2) return 0;
  double mylength = 0;
  for(i=0;i<3;++i) mylength += (ppp[0][i]-ppp[1][i])*(ppp[0][i]-ppp[1][i]);

  return sqrt(mylength);
}
