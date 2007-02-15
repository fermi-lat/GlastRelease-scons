/** @file LinkDef.h
    @brief for rootcint
 $Header$

*/
#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

//#pragma link C++ class VolumeIdentifier+;
#pragma link C++ class McParticle+;
#pragma link C++ class McPositionHit+;
#pragma link C++ class McIntegratingHit+;
#pragma link C++ class McEvent+;
#pragma link C++ class McTrajectory+;
#pragma link C++ class McTrajectoryPoint+;

#endif
