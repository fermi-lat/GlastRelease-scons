#ifndef TAGGERPARAMETERS_HH
#define TAGGERPARAMETERS_HH

// If you change these parameters:
// > make clean
// > make
namespace AncillaryData
{
  const int    N_MODULES              = 4;
  const int    N_MODULES_IN_DATA_WORD = 4;
  const int    N_CHANNELS_PER_MODULE  = 768;
  const int    N_LAYERS_PER_MODULE    = 2;
  const int    N_CHANNELS_PER_LAYER   = 384;
  const int    N_DEBUG_CHANNELS       = 5;
  const int    N_QCD_CHANNELS         = 32;
  const double SIGNAL_THRESHOLD       = 5.0;
  const int    MAX_CLUSTER_GAP        = 5;
  const double STRIPS_PITCH           = 0.0242;
  const double MINIMUM_NOISE          = 1.0;
  const bool   FIT_PEDESTALS          = 0;
}
#endif


