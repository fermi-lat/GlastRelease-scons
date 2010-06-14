/** 
 * @class eventFile::LSE_Info
 *
 * @brief Classes encapsulating event-type-specific acquisition information.
 *
 * @author Bryson Lee <blee@slac.stanford.edu>
 *
 * $Header$
 */

#ifndef LSE_INFO_HH
#define LSE_INFO_HH

#include "eventFile/LSE_GemTime.h"

namespace eventFile {

  struct LSE_Info {
    LSE_Info() {};
    void dump( const char* pre, const char* post ) const;
    unsigned timeTics;
    LSE_GemTime timeHack;
  };

  struct LPA_Info : public LSE_Info {
    LPA_Info() {};
    void dump() const;
    unsigned softwareKey;
    unsigned hardwareKey;
  };

//   struct LCI_Info : public LSE_Info {
//     LCI_Info() {};
//     LCI_Info( const DfiLci::MetaEvent& );
//     void dump( const char* pre, const char* post ) const;
//     bool autoRange;
//     bool zeroSupression;
//     unsigned periodicPrescale;
//     unsigned softwareKey;
//     unsigned writeCfg;
//     unsigned readCfg;
//   };

//   struct LCI_Channel {
//     LCI_Channel() {};
//     LCI_Channel( const DfiLci::Channel& );
//     void dump( const char* pre, const char* post ) const;
//     unsigned short single;
//     bool all;
//     bool latc;
//   };

//   struct LCI_AcdTrigger {
//     LCI_AcdTrigger() {};
//     LCI_AcdTrigger( const DfiLci::AcdTrigger& );
//     void dump( const char* pre, const char* post ) const;
//     unsigned short veto;
//     unsigned short vetoVernier;
//     unsigned short hld;
//   };

//   struct LCI_CalTrigger {
//     LCI_CalTrigger() {};
//     LCI_CalTrigger( const DfiLci::CalTrigger& );
//     void dump( const char* pre, const char* post ) const;
//     unsigned short le;
//     unsigned short he;
//   };

//   struct LCI_ACD_Info : public LCI_Info {
//     ACD_LCI_Info() {};
//     ACD_LCI_Info( const DfiLci::AcdMetaEvent& );
//     void dump( const char* pre, const char* post ) const;
//     unsigned short injected;
//     unsigned short threshold;
//     unsigned short biasDac;
//     unsigned short holdDelay;
//     LCI_AcdTrigger     trigger;
//     LCI_Channel        channel;
//   };

//   struct LCI_CAL_Info : public LCI_Info {
//     CAL_LCI_Info() {};
//     CAL_LCI_Info( const DfiLci::CalMetaEvent& );
//     void dump( const char* pre, const char* post ) const;
//     unsigned short uld;
//     unsigned short injected;
//     unsigned short delay;
//     unsigned short threshold;
//     LCI_CalTrigger     trigger;
//     LCI_Channel        channel;
//   };

//   struct LCI_TKR_Info : public LCI_Info {
//     TKR_LCI_Info() {};
//     TKR_LCI_Info( const DfiLci::TkrMetaEvent& );
//     void dump( const char* pre, const char* post ) const;
//     unsigned short injected;
//     unsigned short delay;
//     unsigned short threshold;
//     Channel        channel;
//   };

};

#endif
