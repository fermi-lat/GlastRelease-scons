/********************************************************************
* mcRootDataCint.h
********************************************************************/
#ifdef __CINT__
#error mcRootDataCint.h/C is only for compilation. Abort cint.
#endif
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define G__ANSIHEADER
#include "G__ci.h"
extern "C" {
extern void G__cpp_setup_tagtablemcRootDataCint();
extern void G__cpp_setup_inheritancemcRootDataCint();
extern void G__cpp_setup_typetablemcRootDataCint();
extern void G__cpp_setup_memvarmcRootDataCint();
extern void G__cpp_setup_globalmcRootDataCint();
extern void G__cpp_setup_memfuncmcRootDataCint();
extern void G__cpp_setup_funcmcRootDataCint();
extern void G__set_cpp_environmentmcRootDataCint();
}


#include "TROOT.h"
#include "TMemberInspector.h"
#include "RMcVertex.h"
#include "RMcParticle.h"
#include "RMcEvent.h"

#ifndef G__MEMFUNCBODY
#endif

extern G__linked_taginfo G__mcRootDataCintLN_TClass;
extern G__linked_taginfo G__mcRootDataCintLN_TObjArray;
extern G__linked_taginfo G__mcRootDataCintLN_TObject;
extern G__linked_taginfo G__mcRootDataCintLN_TLorentzVector;
extern G__linked_taginfo G__mcRootDataCintLN_RMcParticle;
extern G__linked_taginfo G__mcRootDataCintLN_RMcVertex;
extern G__linked_taginfo G__mcRootDataCintLN_RMcEvent;

/* STUB derived class for protected member access */
