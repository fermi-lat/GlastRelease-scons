/** @file LinkDef.h
    @brief for rootcint
 $Header$

*/
#ifdef __CINT__


#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class EventOverlay+;
#pragma link C++ class GlastAxis+;
#pragma link C++ class TkrOverlay+;
#pragma link C++ class AcdOverlay+;
#pragma link C++ class CalOverlay+;
#pragma link C++ class GemOverlayTileList+;
#pragma link C++ class GemOverlayOnePpsTime+;
#pragma link C++ class GemOverlayCondArrivalTime+;
#pragma link C++ class GemOverlay+;
#pragma link C++ class CalDiagDataOverlay+;
#pragma link C++ class TkrDiagDataOverlay+;
#pragma link C++ class DiagDataOverlay+;
#pragma link C++ class PtOverlay+;

#endif

