/** @file LinkDef.h
    @brief for rootcint
 $Header$

*/#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class TVector3+;

#pragma link C++ class ReconEvent+;

#pragma link C++ class AcdRecon+;

#pragma link C++ class CalRecon+;
#pragma link C++ class CalParams+;
#pragma link C++ class CalClusterLayerData+;
#pragma link C++ class CalCluster+;
#pragma link C++ class CalRangeRecData+;
#pragma link C++ class CalXtalRecData+;

#pragma link C++ class TkrRecon+;
#pragma link C++ class TkrCluster+;
#pragma link C++ class TkrTrack+;
#pragma link C++ class TkrTrackHit+;
#pragma link C++ class TkrTrackParams+;
#pragma link C++ class TkrVertex+;

#endif
