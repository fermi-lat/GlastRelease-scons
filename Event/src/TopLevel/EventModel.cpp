// $Id$

#define _GlastEvent_EventModel_CPP_


// Include files
#include "GlastEvent/TopLevel/EventModel.h"
#include "Gaudi/Kernel/Kernel.h"

    
/*! \class EvModel
\brief EvModel class taken from the code supplied by LHCbEvent\TopLevel

  original author Markus Frank, which changes made by Pavel Binko
  this version was originally written by Sawyer Gillespie, hgillesp@u.washington.edu

  Defines the structure of our Transient Data Store
*/
  
class EvModel {
        
public:
    
    EvModel() {
        // The whole GLAST event
        EventModel::Event            = "/Event";
            
        // set up the MC structure
        EventModel::MC::Event = EventModel::Event + "/MC";
        EventModel::MC::McParticle        = EventModel::Event + "/McPaticles";
        EventModel::MC::McIntegratingHits = EventModel::Event + "/McIntegratingHits";
        EventModel::MC::McVertex          = EventModel::Event + "/McVertex";
        EventModel::MC::McPositionHits    = EventModel::Event + "/McPositionHits";

        EventModel::Irf::Event = EventModel::Event + "/Irf";
        EventModel::Irf::IrfAcdHits = EventModel::Irf::Event + "/IrfAcdHits";
        EventModel::Irf::IrfCalHits = EventModel::Irf::Event + "/IrfCalHits";
        EventModel::Irf::IrfTkrHits = EventModel::Irf::Event + "/IrfTkrHits";

        EventModel::Raw::Event = EventModel::Event + "/Raw";
        EventModel::Raw::AcdDigi = EventModel::Raw::Event + "/AcdDigi";
        EventModel::Raw::TdCsIDatas = EventModel::Raw::Event + "/TdCsIDatas";
        EventModel::Raw::TdSiDatas = EventModel::Raw::Event + "/TdSiDatas";

        EventModel::TkrRecon::Event = EventModel::Event + "/TkrRecon";
        EventModel::TkrRecon::SiLayers   = EventModel::TkrRecon::Event + "/SiLayers";
        EventModel::TkrRecon::SiClusters = EventModel::TkrRecon::Event + "/SiClusters";
        EventModel::TkrRecon::SiRecObjs  = EventModel::TkrRecon::Event + "/SiRecObjs";


    }
};
    
    
    static EvModel mod;
    
/*! Class ID definitions for the LHCb Event Model (1st and 2nd byte)
    Reserved numbers are from the GAUDI Framework are 0-99
    Maximum CLID is 65536 = 2^16 - 1

    GLAST POLICY

    MonteCarlo         20xx
    Raw data           21xx
    Digi               22xx
    IrfHit             23xx
    TkrRecon           24xx 
                        -These are classes created by Jose adapted from the 
                         Test Beam Classes.

*/
    
    // we can't redefine the CLID_Event identifier...LHCb already has one of their own
    //const CLID& CLID_Event              = 1000;  // Event root

    //! Irf class IDs
    const CLID& CLID_IrfTkrHit          = 2001;
    const CLID& CLID_IrfAcdHit          = 2002;
    const CLID& CLID_IrfCalHit          = 2003;
    const CLID& CLID_IrfEvent           = 2004;
    const CLID& CLID_IrfTkrLayer        = 2005;


    //! Raw event definitions
    const CLID& CLID_RawEvent           = 2017;
    const CLID& CLID_TdCsIData          = 2014;
    const CLID& CLID_Xtal               = 2015;
    const CLID& CLID_TdCsIDataCnv       = 2016;
    const CLID& CLID_TdSiData           = 2017;
    const CLID& CLID_TdSiDataCnv        = 2018;
    const CLID& CLID_TdGlastData        = 2019;

    //! Monte Carlo class IDs
    const CLID& CLID_MCEvent            = 2012;
    const CLID& CLID_McParticle         = 2100;
    const CLID& CLID_McIntegratingHit   = 2101;
    const CLID& CLID_McVertex           = 2102;
    const CLID& CLID_McPositionHit      = 2103;

    //! Digitization class IDs
    const CLID& CLID_AcdDigi            = 2201;

    
    //! Classes adapted by Jose from the test beam
    const CLID& CLID_SiLayers           = 2452;
    const CLID& CLID_SiClusters         = 2454;
    const CLID& CLID_SiRecObjs          = 2456;
        
