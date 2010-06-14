#ifndef TAGGERMODULE_HH
#define TAGGERMODULE_HH

#include "AncillaryDataEvent/TaggerLayer.h"

class TaggerModule

{
 public:
 
  TaggerModule(int moduleId);
  ~TaggerModule();

  void processData();
  void resetData();
  
  int getId()
    const {return m_moduleId;}
  TaggerLayer *getLayer(int layerId)
    const {return m_layers[layerId];}
  TaggerChannel *getChannel(int layerId, int channelId)
    const {return getLayer(layerId)->getChannel(channelId);}
  double getXCrossingPosition() const {return m_xCrossingPosition;}
  double getXCrossingError()    const {return m_xCrossingError;}
  double getYCrossingPosition() const {return m_yCrossingPosition;}
  double getYCrossingError()    const {return m_yCrossingError;}
  
 private:
  
  int  m_moduleId;
  TaggerLayer *m_layers[N_LAYERS_PER_MODULE];
  double m_xCrossingPosition;
  double m_xCrossingError;
  double m_yCrossingPosition;
  double m_yCrossingError;
};
#endif
