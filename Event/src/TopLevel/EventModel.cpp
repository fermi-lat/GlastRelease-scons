// File and Version Information:
// $Header$

#define _Event_EventModel_CPP_


#include "Event/TopLevel/EventModel.h"
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/ClassID.h"

    
/** @class EvModel
 *  @brief Event Model: Definition of logical paths and class identifiers
 *
 * $Header$
 */
class EvModel {
        
public:
    
    EvModel() {
        // Access to GLAST event
        EventModel::EventHeader              = "/Event";
            
        // Monte Carlo 
        EventModel::MC::Event                = EventModel::EventHeader + "/MC";
        EventModel::MC::McParticleCol        = EventModel::MC::Event  + "/McParticleCol";

        EventModel::MC::McPositionHitCol     = EventModel::MC::Event  + "/PositionHitsCol";
        EventModel::MC::McIntegratingHitCol  = EventModel::MC::Event  + "/IntegratingHitsCol";
        EventModel::MC::McTkrStripCol        = EventModel::MC::Event  + "/StripCol";
        EventModel::MC::D2EntryCol           = EventModel::MC::Event  + "/D2EntryCol";
        EventModel::MC::ExposureCol          = EventModel::MC::Event  + "/ExposureCol";

        EventModel::MC::McEventStructure     = EventModel::MC::Event  + "/McEventStructure";

        EventModel::MC::McPartToHitTab       = EventModel::MC::Event  + "/McPartToHitTab";
        EventModel::MC::McClusToLyrHitTab    = EventModel::MC::Event  + "/McClusToLyrHitTab";
        EventModel::MC::McLyrToHitTab        = EventModel::MC::Event  + "/McLyrToHitTab";
        EventModel::MC::McSiLayerHitCol      = EventModel::MC::Event  + "/McSiLayerHitCol";

        // Digi event
        EventModel::Digi::Event              = EventModel::EventHeader + "/Digi";
        EventModel::Digi::AcdDigiCol         = EventModel::Digi::Event + "/AcdDigiCol";
        EventModel::Digi::TkrDigiCol         = EventModel::Digi::Event + "/TkrDigiCol";
        EventModel::Digi::CalDigiCol         = EventModel::Digi::Event + "/CalDigiCol";
        EventModel::Digi::CalDigiHitTab      = EventModel::Digi::Event + "/CalDigiHitTab";
        EventModel::Digi::TkrDigiHitTab      = EventModel::Digi::Event + "/TkrDigiHitTab";
        EventModel::Digi::TkrClusterHitTab   = EventModel::Digi::Event + "/TkrClusterHitTab";
        

        // reconstructed data (Tracker)
        EventModel::TkrRecon::Event          = EventModel::EventHeader + "/TkrRecon";
        EventModel::TkrRecon::SiLayers       = EventModel::TkrRecon::Event + "/SiLayers";
        EventModel::TkrRecon::TkrClusterCol  = EventModel::TkrRecon::Event + "/TkrClusterCol";
        EventModel::TkrRecon::TkrPatCandCol  = EventModel::TkrRecon::Event + "/TkrPatCandCol";
        EventModel::TkrRecon::SiRecObjs      = EventModel::TkrRecon::Event + "/SiRecObjs";
        EventModel::TkrRecon::TkrFitTrackCol = EventModel::TkrRecon::Event + "/TkrFitTrackCol";
        EventModel::TkrRecon::TkrTrackTab    = EventModel::TkrRecon::Event + "/TkrTrackTab";
        EventModel::TkrRecon::TkrVertexCol   = EventModel::TkrRecon::Event + "/TkrVertexCol";
        EventModel::TkrRecon::TkrVertexTab   = EventModel::TkrRecon::Event + "/TkrVertexTab";


        //reconstructed Cal data

		EventModel::CalRecon::Event          = EventModel::EventHeader + "/CalRecon";
		EventModel::CalRecon::CalXtalRecCol  = EventModel::CalRecon::Event + "/CalXtalRecCol";
		EventModel::CalRecon::CalClusterCol  = EventModel::CalRecon::Event + "/CalClusterCol";

		
		
		// reconstructed ACD data
        EventModel::AcdRecon::Event          = EventModel::EventHeader + "/AcdRecon";
    }
};
    
    
    static EvModel mod;    // where  used? has file scope     
    
/*  Class ID definitions for the Glast Event Model
    Maximum CLID is 65536 = 2^16 - 1

    Category           ID range      Comment
    ---------          --------      -------
    Gaudi Kernel          0 -  99    see Gaudi/Kernel/Kernel.cpp
                        100 - 199    general classes (Run, Event, EventTag, ContainedObject)
                                     see Gaudi/Kernel/Kernel.cpp
                                     Some of these classes may be removed in the next 
                                     Gaudi release

    GLAST class Categories and ID ranges (proposal)
    
    Category           ID range      Comment
    ---------          --------      -------
    EventSelection     200  -  299   High level event information used for event selection
                                     (Run and EventTag are at present in Kernel.cpp)
    Utilities          300  -  399   used at several places of the EventModel
    
    MonteCarlo         1000
      Generator        1001 - 1099   Flux generator related classes
      Kine+Hit         1100 - 1199   Kinematics, Hits and associations
      IrfHit           1200 - 1299   Special IRF Hits

    Digi data          1300 - 1399   Digi (and possible Hit/Digi associations)

    Raw data           1400 - 1499   Raw data (and possible Hit/Digi associations)

    Trigger (simulation)
      Trigger          2000 - 2099   Trigger summary 
      LVL1             2100 - 2199   LVL1 classes
      Hlt              2200 - 2299   Higher level trigger classes

    Reconstruction
      Recon            3000 - 3099   Combined Recon and Recon summary
      TkrRecon         3100 - 3199   Tracker recon
      CalRecon         3200 - 3299   Calorimeter recon
      AcdRecon         3300 - 3399   ACD recon

    Analysis           4000 - 4999   Event interpretation beyond reconstruction, e.g.
                                     information for event classification, diagnosis, ... 
                                     Typically not written to persistent store.
    Note:
    Detector description and calibrations have their own trees, not connected
    to the EventModel tree.   Calibration CLID's will be contained within
    the range 6000 - 6999.
 */
    
    // Declaration of Identifiers
    // The order is: General, Tracker, Calorimeter, ACD
 
    // HMK No longer defined in Gaudi - it is left to the user to define
    const CLID& CLID_Event            =  110;

    //! Monte Carlo class IDs
    const CLID& CLID_McEvent            = 1100;
    const CLID& CLID_McVertex           = 1101;
    const CLID& CLID_McParticle         = 1102;
    const CLID& CLID_McParticleCol      = CLID_McParticle+CLID_ObjectList;
    const CLID& CLID_McPositionHit      = 1103;
    const CLID& CLID_McPositionHitCol   = CLID_McPositionHit+CLID_ObjectVector;
    const CLID& CLID_McIntegratingHit   = 1104;
    const CLID& CLID_McIntegratingHitCol= CLID_McIntegratingHit+CLID_ObjectVector;
    const CLID& CLID_McTrajectory       = 1105;
    const CLID& CLID_McTkrStrip         = 1106;
    const CLID& CLID_D2Entry            = 1107;
    const CLID& CLID_Exposure           = 1108;

    //! Raw event and Digi IDs
    const CLID& CLID_DigiEvent          = 1300; 
    const CLID& CLID_TkrDigi            = 1310;  // indicative only, use 1310, 1311,.. for Tkr
    const CLID& CLID_TkrDigiCol         = CLID_TkrDigi+CLID_ObjectVector;
    const CLID& CLID_CalDigi            = 1320;  // indicative only, use 1320, 1321,.. for Cal
    const CLID& CLID_CalDigiCol         = CLID_CalDigi+CLID_ObjectVector;
    const CLID& CLID_AcdDigi            = 1330;
    const CLID& CLID_AcdDigiCol         = CLID_AcdDigi+CLID_ObjectVector;
    
    // clarify where these are used 
    const CLID& CLID_RawEvent           = 1400;  // temporary...will be replaced by DataEvent
    const CLID& CLID_TdGlastData        = 1401;
    const CLID& CLID_TdSiData           = 1411;
    const CLID& CLID_TdSiDataCnv        = 1412;
    const CLID& CLID_TdCsIData          = 1421;
    const CLID& CLID_TdCsIDataCnv       = 1422;
    const CLID& CLID_Xtal               = 1423;
    const CLID& CLID_TdVetoData         = 1424;
    const CLID& CLID_LdGlastData        = 1425;
    const CLID& CLID_LdGlastDataCnv     = 1426;
   
    //! Reconstruction
    const CLID& CLID_RecEvent           = 3000;
  
    //! Reconstruction: Tkr class IDs
    const CLID& CLID_TkrRecon           = 3100;
    const CLID& CLID_SiLayers           = 3101;  // Use Tkr instead of Si ?
    const CLID& CLID_TkrClusterCol      = 3102;  // Is this really Layers, Clusters, RecObjs
    const CLID& CLID_TkrPatCandCol      = 3103;  // Is this really Layers, Clusters, RecObjs
    const CLID& CLID_SiRecObjs          = 3104;  // or rather Layer, Cluster, RecObj objects?
    const CLID& CLID_TkrFitTrackCol     = 3105;  // or rather Layer, Cluster, RecObj objects?
    const CLID& CLID_TkrVertex          = 3106;  
    const CLID& CLID_TkrVertexCol       = CLID_TkrVertex+CLID_ObjectVector;// CLID_TkrVertexCol is ObjectVector of TkrVertices

    //! Reconstruction: Cal class IDs
    const CLID& CLID_CalRecon           = 3200;
    const CLID& CLID_CalXtalRecData     = 3201;
    const CLID& CLID_CalXtalRecCol      = CLID_CalXtalRecData+CLID_ObjectVector;
    const CLID& CLID_CalClusterCol      = 3202;

    //! Reconstruction: Acd class IDs
    const CLID& CLID_AcdRecon           = 3300;
        

    //! Analysis
    const CLID& CLID_AnalEvent          = 4000;
        
    //! Classes adapted by Sasha Chekhtman from tb calorimeter reconstruction    
    const CLID& CLID_CalADCLogs         = 2601;
    const CLID& CLID_CalRecLogs         = 2602;

    //! Utilities
    const CLID& CLID_RefTable1to1       =  321;
    const CLID& CLID_RefTable1toN       =  322;

