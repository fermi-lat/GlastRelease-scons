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
#include "AcdDigi.h"
#include "AcdTile.h"
#include "CalDigi.h"
#include "CalHit.h"
#include "ESAPID.h"
#include "Event.h"
#include "L1T.h"
#include "LogID.h"
#include "StripID.h"
#include "Tagger.h"
#include "TileID.h"
#include "TkrDigi.h"
#include "TkrLayer.h"

#ifndef G__MEMFUNCBODY
#endif

extern G__linked_taginfo G__EventCintLN_TClass;
extern G__linked_taginfo G__EventCintLN_TObjArray;
extern G__linked_taginfo G__EventCintLN_TObject;
extern G__linked_taginfo G__EventCintLN_TClonesArray;
extern G__linked_taginfo G__EventCintLN_TileID;
extern G__linked_taginfo G__EventCintLN_AcdTile;
extern G__linked_taginfo G__EventCintLN_AcdTilecLcLdA;
extern G__linked_taginfo G__EventCintLN_AcdDigi;
extern G__linked_taginfo G__EventCintLN_LogID;
extern G__linked_taginfo G__EventCintLN_LogIDcLcLdA;
extern G__linked_taginfo G__EventCintLN_LogIDcLcLCALAxes;
extern G__linked_taginfo G__EventCintLN_CalHit;
extern G__linked_taginfo G__EventCintLN_CalHitcLcLdA;
extern G__linked_taginfo G__EventCintLN_CalHitcLcLADCGain;
extern G__linked_taginfo G__EventCintLN_CalHitcLcLLogEnd;
extern G__linked_taginfo G__EventCintLN_CalDigi;
extern G__linked_taginfo G__EventCintLN_ESAPID;
extern G__linked_taginfo G__EventCintLN_L1T;
extern G__linked_taginfo G__EventCintLN_StripID;
extern G__linked_taginfo G__EventCintLN_StripIDcLcLdA;
extern G__linked_taginfo G__EventCintLN_TkrLayer;
extern G__linked_taginfo G__EventCintLN_TkrLayercLcLTKRAxes;
extern G__linked_taginfo G__EventCintLN_TkrDigi;
extern G__linked_taginfo G__EventCintLN_Event;
extern G__linked_taginfo G__EventCintLN_Tagger;

/* STUB derived class for protected member access */
