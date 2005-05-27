/** @file LinkDef.h
    @brief for rootcint
 $Header$

*/
#ifdef __CINT__


#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class Event+;
#pragma link C++ class EventDisplay+;
#pragma link C++ class Layer+;
#pragma link C++ class Tracker+;
#pragma link C++ class Recon+;
#pragma link C++ class Efficiency+; 
#pragma link C++ class TriggerEfficiency+; 
#pragma link C++ class Residual+;

#endif

