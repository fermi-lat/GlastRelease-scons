// $Id$

#define _GlastEvent_EventModel_CPP_


// Include files
#include "GlastEvent/TopLevel/EventModel.h"
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/ClassID.h"

    
/*! \class EvModel
 *  \brief Event Model: Definition of logical paths and class identifiers
 *
 */
  
class EvModel {
        
public:
    
    EvModel() {
        // Access to GLAST event
        EventModel::Event                 = "/Event";
            
        // Monte Carlo 
        EventModel::MC::Event             = EventModel::Event + "/MC";
        EventModel::MC::McVertexCol        = EventModel::MC::Event  + "/McVertexCol";
        EventModel::MC::McParticleCol       = EventModel::MC::Event  + "/McParticleCol";

        EventModel::MC::McPositionHits    = EventModel::MC::Event  + "/McPositionHits";
        EventModel::MC::McIntegratingHits = EventModel::MC::Event  + "/McIntegratingHits";

        // Digi event
        EventModel::Digi::Event           = EventModel::Event + "/Digi";
        EventModel::Digi::AcdDigis        = EventModel::Digi::Event + "/AcdDigis";
        EventModel::Digi::TkrDigis        = EventModel::Digi::Event + "/TkrDigis";

        // Irf event
        EventModel::Irf::Event            = EventModel::Event + "/Irf";
        EventModel::Irf::IrfTkrHits       = EventModel::Irf::Event + "/IrfTkrHits";
        EventModel::Irf::IrfCalHits       = EventModel::Irf::Event + "/IrfCalHits";
        EventModel::Irf::IrfAcdHits       = EventModel::Irf::Event + "/IrfAcdHits";

        // Data Data
        EventModel::Data::Event           = EventModel::Event + "/Data";
        EventModel::Data::TdGlastData     = EventModel::Data::Event + "/TdGlastData";
        EventModel::Data::TdSiData        = EventModel::Data::Event + "/TdSiData";
        EventModel::Data::TdCsIData       = EventModel::Data::Event + "/TdCsIData";
        EventModel::Data::TdVetoData      = EventModel::Data::Event + "/TdVetoData";
        

        // reconstructed data (Tracker)
        EventModel::TkrRecon::Event       = EventModel::Event + "/TkrRecon";
        EventModel::TkrRecon::SiLayers    = EventModel::TkrRecon::Event + "/SiLayers";
        EventModel::TkrRecon::SiClusters  = EventModel::TkrRecon::Event + "/SiClusters";
        EventModel::TkrRecon::SiRecObjs   = EventModel::TkrRecon::Event + "/SiRecObjs";

        // reconstructed ACD data
        EventModel::AcdRecon::Event       = EventModel::Event + "/AcdRecon";
    }
};
    
    
    static EvModel mod;    // where  used? has file scope     
    
/*  Class ID definitions for the Glast Event Model
    Maximum CLID is 65536 = 2^16 - 1

    Categorie          ID range      Comment
    ---------          --------      -------
    Gaudi Kernel          0 -  99    see Gaudi/Kernel/Kernel.cpp
                        100 - 199    general classes (Run, Event, EventTag, ContainedObject)
                                     see Gaudi/Kernel/Kernel.cpp
                                     Some of these classes may be removed in the next 
                                     Gaudi release

    GLAST class Categories and ID ranges (proposal)
    
    Categorie          ID range      Comment
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
    to the EventModel tree.  
 */
    
    // Declaration of Identifiers
    // The order is: General, Tracker, Calorimeter, ACD
 
    //const CLID& CLID_Event            =  110;  // defined in Gaudi/Kernel/Kernel.cpp

    //! Monte Carlo class IDs
    const CLID& CLID_McEvent            = 1100;
    const CLID& CLID_McVertex           = 1101;
    const CLID& CLID_McParticle         = 1102;
    const CLID& CLID_McPositionHit      = 1103;
    const CLID& CLID_McIntegratingHit   = 1104;
    const CLID& CLID_McVertexCnv        = 1105;

    //! Irf class IDs
    const CLID& CLID_IrfEvent           = 1200;
    const CLID& CLID_IrfTkrHit          = 1210;
    const CLID& CLID_IrfTkrLayer        = 1211;
    const CLID& CLID_IrfCalHit          = 1220;
    const CLID& CLID_IrfAcdHit          = 1230;

    //! Raw event and Digi IDs
    const CLID& CLID_DigiEvent          = 1300; 
    const CLID& CLID_TkrDigi            = 1310;  // indicative only, use 1310, 1311,.. for Tkr
    const CLID& CLID_CalDigi            = 1320;  // indicative only, use 1320, 1321,.. for Cal
    const CLID& CLID_AcdDigi            = 1330;
    
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
    const CLID& CLID_SiClusters         = 3102;  // Is this really Layers, Clusters, RecObjs
    const CLID& CLID_SiRecObjs          = 3103;  // or rather Layer, Cluster, RecObj objects?

    //! Reconstruction: Tkr class IDs
    const CLID& CLID_CalRecon           = 3200;

    //! Reconstruction: Acd class IDs
    const CLID& CLID_AcdRecon           = 3300;
        

    //! Analysis
    const CLID& CLID_AnalEvent          = 4000;
        
    //! Classes adapted by Sasha Chekhtman from tb calorimeter reconstruction
    
    const CLID& CLID_CalADCLogs         = 2601;
    const CLID& CLID_CalRecLogs         = 2602;
    const CLID& CLID_CalClusterList     = 2603;

    //! Utilities
    const CLID& CLID_RefTable1to1       =  321;
    const CLID& CLID_RefTable1toN       =  322;

