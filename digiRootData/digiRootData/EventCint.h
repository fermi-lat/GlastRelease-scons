/********************************************************************
* EventCint.h
********************************************************************/
#ifdef __CINT__
#error EventCint.h/C is only for compilation. Abort cint.
#endif
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define G__ANSIHEADER
#include "G__ci.h"
extern "C" {
extern void G__cpp_setup_tagtableEventCint();
extern void G__cpp_setup_inheritanceEventCint();
extern void G__cpp_setup_typetableEventCint();
extern void G__cpp_setup_memvarEventCint();
extern void G__cpp_setup_globalEventCint();
extern void G__cpp_setup_memfuncEventCint();
extern void G__cpp_setup_funcEventCint();
extern void G__set_cpp_environmentEventCint();
}


#include "TROOT.h"
#include "TMemberInspector.h"
#include "AcdHeader.h"
#include "AcdId.h"
#include "AcdDigi.h"
#include "CalHeader.h"
#include "CalDigi.h"
#include "CalLogId.h"
#include "ESAPID.h"
#include "Event.h"
#include "L1T.h"
#include "StripId.h"
#include "SubSystemHeader.h"
#include "Tagger.h"
#include "TkrHeader.h"
#include "TkrDigi.h"

#ifndef G__MEMFUNCBODY
#endif

extern G__linked_taginfo G__EventCintLN_bool;
extern G__linked_taginfo G__EventCintLN_TClass;
extern G__linked_taginfo G__EventCintLN_TObjArray;
extern G__linked_taginfo G__EventCintLN_TObject;
extern G__linked_taginfo G__EventCintLN_SubSystemHeader;
extern G__linked_taginfo G__EventCintLN_AcdHeader;
extern G__linked_taginfo G__EventCintLN_AcdId;
extern G__linked_taginfo G__EventCintLN_AcdIdcLcLdA;
extern G__linked_taginfo G__EventCintLN_AcdDigi;
extern G__linked_taginfo G__EventCintLN_AcdDigicLcLdA;
extern G__linked_taginfo G__EventCintLN_CalHeader;
extern G__linked_taginfo G__EventCintLN_CalLogId;
extern G__linked_taginfo G__EventCintLN_CalLogIdcLcLdA;
extern G__linked_taginfo G__EventCintLN_CalLogIdcLcLCALAxes;
extern G__linked_taginfo G__EventCintLN_CalDigi;
extern G__linked_taginfo G__EventCintLN_CalDigicLcLAdcRange;
extern G__linked_taginfo G__EventCintLN_CalDigicLcLLogFace;
extern G__linked_taginfo G__EventCintLN_CalDigicLcLCalTrigMode;
extern G__linked_taginfo G__EventCintLN_CalDigicLcLdA;
extern G__linked_taginfo G__EventCintLN_ESAPID;
extern G__linked_taginfo G__EventCintLN_L1T;
extern G__linked_taginfo G__EventCintLN_TClonesArray;
extern G__linked_taginfo G__EventCintLN_TkrHeader;
extern G__linked_taginfo G__EventCintLN_StripId;
extern G__linked_taginfo G__EventCintLN_StripIdcLcLdA;
extern G__linked_taginfo G__EventCintLN_TkrDigi;
extern G__linked_taginfo G__EventCintLN_TkrDigicLcLTKRAxes;
extern G__linked_taginfo G__EventCintLN_Event;
extern G__linked_taginfo G__EventCintLN_Tagger;

/* STUB derived class for protected member access */
