/** @file LinkDef.h
    @brief for rootcint
 $Header$

*/
#ifdef __CINT__


#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class DigiEvent+;
#pragma link C++ class GlastAxis+;
#pragma link C++ class TkrDigi+;
#pragma link C++ class L1T+;
#pragma link C++ class AcdDigi+;
#pragma link C++ class CalDigi+;
#pragma link C++ class CalXtalReadout+;
#pragma link C++ class TkrDiagnosticData+; 
#pragma link C++ class CalDiagnosticData+;

/// LSF stuff
#pragma link C++ class GemTime+;  
#pragma link C++ class TimeTone+;  
#pragma link C++ class LsfTime+;  
#pragma link C++ class Configuration+;  
#pragma link C++ class GemScalers+; 
#pragma link C++ class RunInfo+;  
#pragma link C++ class DatagramInfo+; 
#pragma link C++ class Configuration+;  
#pragma link C++ class LpaConfiguration+; 
#pragma link C++ class Channel+; 
#pragma link C++ class LciConfiguration+; 
#pragma link C++ class LciAcdConfiguration+;
#pragma link C++ class LciCalConfiguration+;
#pragma link C++ class LciTkrConfiguration+;
#pragma link C++ class LciAcdConfiguration::AcdTrigger+;
#pragma link C++ class LciCalConfiguration::CalTrigger+;
#pragma link C++ class MetaEvent+;
#pragma link C++ class Ccsds+;

#pragma link C++ class EventSummaryData+;
#pragma link C++ class GemTileList+;
#pragma link C++ class GemOnePpsTime+;
#pragma link C++ class GemCondArrivalTime+;
#pragma link C++ class Gem+;
#pragma link C++ class ErrorData+;
#pragma link C++ class Tem+;

// OBF
#pragma link C++ class CalLogInfo+;
#pragma link C++ class TfcProjection+;
#pragma link C++ class TfcProjectionDir+;
#pragma link C++ class TfcProjectionCol+;
#pragma link C++ class ObfTrack+;
#pragma link C++ class FilterStatus+;

// Beamtest Ancillary Data 2006
#pragma link C++ class TaggerHit+;
#pragma link C++ class QdcHit+;
#pragma link C++ class AdfDigi+;

#endif

