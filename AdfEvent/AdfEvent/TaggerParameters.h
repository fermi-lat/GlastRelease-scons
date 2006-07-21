#ifndef TAGGERPARAMETERS_HH
#define TAGGERPARAMETERS_HH

// If you change these parameters:
// > make clean
// > make
namespace AncillaryData
{
  // ALL THE UNITS ARE IN mm
  const unsigned int    N_MODULES              = 4;
  const unsigned int    N_MODULES_IN_DATA_WORD = 4;
  const unsigned int    N_CHANNELS_PER_MODULE  = 768;
  const unsigned int    N_LAYERS_PER_MODULE    = 2;
  const unsigned int    N_CHANNELS_PER_LAYER   = 384;
  const unsigned int    N_DEBUG_CHANNELS       = 5;
  const unsigned int    N_QDC_MODULES          = 1;
  const unsigned int    N_QDC_CHANNELS         = 32;
  const unsigned int  MAX_CLUSTER_GAP          =5;
  const double SIGNAL_THRESHOLD       = 5.0;
  const double STRIPS_PITCH           = 0.242;
  const double LAYER_WIDTH            = 0.400;
  const double MINIMUM_NOISE          = 1.0;
  const bool   FIT_PEDESTALS          = 1;
}
#endif


