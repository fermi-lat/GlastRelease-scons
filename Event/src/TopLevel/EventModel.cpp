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

        EventModel::Irf::Event = EventModel::Event + "/Irf";
        EventModel::Irf::IrfAcdHits = EventModel::Irf::Event + "/IrfAcdHits";
        EventModel::Irf::IrfCalHits = EventModel::Irf::Event + "/IrfCalHits";
        EventModel::Irf::IrfTkrHits = EventModel::Irf::Event + "/IrfTkrHits";

        EventModel::Raw::Event = EventModel::Event + "/Raw";
        EventModel::Raw::TdCsIDatas = EventModel::Raw::Event + "/TdCsIDatas";
        EventModel::Raw::TdSiDatas = EventModel::Raw::Event + "/TdSiDatas";
    }
};
    
    
    static EvModel mod;
    
/*! Class ID definitions for the LHCb Event Model (1st and 2nd byte)
    Reserved numbers are from the GAUDI Framework are 0-99
    Maximum CLID is 65536 = 2^16 - 1
*/
    
    // we can't redefine the CLID_Event identifier...LHCb already has one of their own
    //const CLID& CLID_Event              = 1000;  // Event root

    // Irf class IDs
    const CLID& CLID_IrfTkrHit          = 2001;
    const CLID& CLID_IrfAcdHit          = 2002;
    const CLID& CLID_IrfCalHit          = 2003;
    const CLID& CLID_IrfEvent           = 2004;
    const CLID& CLID_IrfTkrLayer        = 2005;

    // Monte Carlo class IDs
  //  const CLID& CLID_MCTKRHit           = 2008;  
  //  const CLID& CLID_MCACDHit           = 2009;
    const CLID& CLID_MCTrack            = 2010; 
 //   const CLID& CLID_MCCalorimeterHit   = 2011;
    const CLID& CLID_MCEvent            = 2012;
//    const CLID& CLID_MCSiLayer          = 2013;

    //! Raw event definitions
    const CLID& CLID_RawEvent           = 2017;
    const CLID& CLID_TdCsIData            = 2014;
    const CLID& CLID_Xtal               = 2015;
    const CLID& CLID_TdCsIDataCnv         = 2016;
    const CLID& CLID_TdSiData           =  2017;
    const CLID& CLID_TdSiDataCnv        =   2018;
    const CLID& CLID_TdGlastData        =   2019;
        
