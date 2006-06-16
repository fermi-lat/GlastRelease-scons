#ifndef TAGGERCHANNEL_HH
#define TAGGERCHANNEL_HH

// Standard stuff.
#include <stdlib.h>
#include <iostream>
#include <fstream>

// Our own classes.
#include "TaggerParameters.h"

class TaggerChannel

{
 public:
 
  TaggerChannel(int channelId, int layerId, int moduleId);
  ~TaggerChannel(){;}
  
  void processData();
  void reset();
  void subtractPedestal();
  
  int    getId()            const {return m_channelId;}
  double getPedestal()      const {return m_pedestal;}
  double getRawNoise()      const {return m_rawNoise;}
  double getCmsNoise()      const {return m_cmsNoise;}
  bool   isBad()            const {return m_isBad;}
  int    getRawADC()        const {return m_rawADC;}
  double getSubtractedADC() const {return m_subtractedADC;}
  bool   isHit()            const {return m_isHit;}
  
  void setPedestal(float pedestal)  {m_pedestal = pedestal;}
  void setRawNoise(float raw_noise) {m_rawNoise = raw_noise;}
  void setCmsNoise(float cms_noise) {m_cmsNoise = cms_noise;}
  void setBadFlag(bool badFlag)     {m_isBad = badFlag;}
  void setRawADC(int adc)           {m_rawADC = adc;}
 
 private:
 
  bool   m_DEBUG;
  int    m_channelId;
  int    m_layerId;
  int    m_moduleId;
  double m_pedestal;
  double m_rawNoise;
  double m_cmsNoise;
  bool   m_isBad;
  int    m_rawADC;
  double m_subtractedADC;
  bool   m_isHit;
};
#endif
