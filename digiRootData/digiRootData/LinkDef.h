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
//#pragma link C++ class TowerId+;
#pragma link C++ class TkrDigi+;
#pragma link C++ class L1T+;
//#pragma link C++ class AcdId+;
#pragma link C++ class AcdDigi+;
#pragma link C++ class CalDigi+;
#pragma link C++ class CalXtalReadout+;
//#pragma link C++ class CalXtalId+;
#pragma link C++ class TkrDiagnosticData+; //class name different from filename!
#pragma link C++ class CalDiagnosticData+; //class name different from filename!
#pragma link C++ class EventSummaryData+;
#pragma link C++ class GemTileList+;
#pragma link C++ class GemOnePpsTime+;
#pragma link C++ class GemCondArrivalTime+;
#pragma link C++ class Gem+;
#pragma link C++ class ErrorData+;
#pragma link C++ class Tem+;

#endif

