/** -*- Mode: C++ -*-
 * @class eventFile::LSE_Info
 *
 * @brief Classes encapsulating event-type-specific acquisition information.
 *
 * @author Bryson Lee <blee@slac.stanford.edu>
 *
 * $Header$
 */

#ifndef EVENTFILE_LSE_INFO_HH
#define EVENTFILE_LSE_INFO_HH

#include <vector>

#include "eventFile/LSE_GemTime.h"

namespace eventFile {

  class LPA_Handler;
  class LSEReader;
  class LSEWriter;

  struct LSE_Info {
    typedef enum _InfoType {
      NONE = -1,
      LPA  =  0,
      LCI_ACD,
      LCI_CAL,
      LCI_TKR,
      NumLSEInfoTypes,
    } InfoType;

    LSE_Info() {};
    void dump( const char* pre, const char* post ) const;
    unsigned timeTics;
    LSE_GemTime timeHack;
    int compressionLevel;
    int compressedSize;
  };

  struct LPA_Info : public LSE_Info {
    LPA_Info() {};
    void dump() const;
    unsigned softwareKey;
    unsigned hardwareKey;
    unsigned lpaDbKey;
    std::vector< LPA_Handler > handlers;

  private:
    void write( FILE* fp ) const;
    void read( FILE* fp );
    friend class LSEReader;
    friend class LSEWriter;

  };

  struct LCI_Info : public LSE_Info {
    LCI_Info() {};
    void dump( const char* pre, const char* post ) const;
    bool autoRange;
    bool zeroSupression;
    bool strobe;
    unsigned periodicPrescale;
    unsigned softwareKey;
    unsigned writeCfg;
    unsigned readCfg;
  };

  struct LCI_Channel {
    LCI_Channel() {};
    LCI_Channel( unsigned short n, bool s, bool a, bool l, bool p ) :
      numChan(n), single(s), all(a), latc(l), perFE(p) {};
    void dump( const char* pre, const char* post ) const;
    unsigned short numChan;
    bool single;
    bool all;
    bool latc;
    bool perFE;
  };

  struct LCI_AcdTrigger {
    LCI_AcdTrigger() {};
    LCI_AcdTrigger( unsigned short v, unsigned short vV, unsigned short h ) :
      veto(v), vetoVernier(vV), hld(h) {};
    void dump( const char* pre, const char* post ) const;
    unsigned short veto;
    unsigned short vetoVernier;
    unsigned short hld;
  };

  struct LCI_CalTrigger {
    LCI_CalTrigger() {};
    LCI_CalTrigger( unsigned short l, unsigned short lTE, unsigned short h, unsigned short hTE ) :
      le(l), lowTrgEna(lTE), he(h),  highTrgEna(hTE) {};
    void dump( const char* pre, const char* post ) const;
    unsigned short le;
    unsigned short lowTrgEna;
    unsigned short he;
    unsigned short highTrgEna;
  };

  struct LCI_ACD_Info : public LCI_Info {
    LCI_ACD_Info() {};
    void dump() const;
    unsigned short injected;
    unsigned short threshold;
    unsigned short biasDac;
    unsigned short holdDelay;
    unsigned short hitmapDelay;
    unsigned short range;
    LCI_AcdTrigger     trigger;
    LCI_Channel        channel;
  };

  struct LCI_CAL_Info : public LCI_Info {
    LCI_CAL_Info() {};
    void dump() const;
    unsigned short uld;
    unsigned short injected;
    unsigned short delay;
    unsigned short firstRange;
    unsigned short threshold;
    unsigned short calibGain;
    unsigned short highCalEna;
    unsigned short highRngEna;
    unsigned short highGain;
    unsigned short lowCalEna;
    unsigned short lowRngEna;
    unsigned short lowGain;
    LCI_CalTrigger     trigger;
    LCI_Channel        channel;
  };

  struct LCI_TKR_Info : public LCI_Info {
    LCI_TKR_Info() {};
    void dump() const;
    enum DacRange { DAC_LO=0, DAC_HI=1 };
    DacRange dacRange;
    unsigned short injected;
    unsigned short delay;
    unsigned short threshold;
    unsigned short splitLow;
    unsigned short splitHigh;
    LCI_Channel        channel;
  };

};

#endif
