/********************************************************************
* ReconCint.h
********************************************************************/
#ifdef __CINT__
#error ReconCint.h/C is only for compilation. Abort cint.
#endif
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define G__ANSIHEADER
#include "G__ci.h"
extern "C" {
extern void G__cpp_setup_tagtableReconCint();
extern void G__cpp_setup_inheritanceReconCint();
extern void G__cpp_setup_typetableReconCint();
extern void G__cpp_setup_memvarReconCint();
extern void G__cpp_setup_globalReconCint();
extern void G__cpp_setup_memfuncReconCint();
extern void G__cpp_setup_funcReconCint();
extern void G__set_cpp_environmentReconCint();
}


#include "TROOT.h"
#include "TMemberInspector.h"
#include "CalCluster.h"
#include "CalLogEne.h"
#include "CalRecon.h"
#include "Recon.h"
#include "ReconHeader.h"
#include "TkrGamma.h"
#include "TkrHit.h"
#include "TkrLocator.h"
#include "TkrRecon.h"
#include "TkrSiCluster.h"
#include "TkrTrack.h"

#ifndef G__MEMFUNCBODY
#endif

extern G__linked_taginfo G__ReconCintLN_TClass;
extern G__linked_taginfo G__ReconCintLN_TObjArray;
extern G__linked_taginfo G__ReconCintLN_TObject;
extern G__linked_taginfo G__ReconCintLN_TVector2;
extern G__linked_taginfo G__ReconCintLN_TVector3;
extern G__linked_taginfo G__ReconCintLN_TArrayF;
extern G__linked_taginfo G__ReconCintLN_CalCluster;
extern G__linked_taginfo G__ReconCintLN_LogId;
extern G__linked_taginfo G__ReconCintLN_CalLogEne;
extern G__linked_taginfo G__ReconCintLN_CalRecon;
extern G__linked_taginfo G__ReconCintLN_TkrGamma;
extern G__linked_taginfo G__ReconCintLN_TkrTrack;
extern G__linked_taginfo G__ReconCintLN_TkrSiCluster;
extern G__linked_taginfo G__ReconCintLN_TkrRecon;
extern G__linked_taginfo G__ReconCintLN_ReconHeader;
extern G__linked_taginfo G__ReconCintLN_ReconHeadercLcLdA;
extern G__linked_taginfo G__ReconCintLN_ReconHeadercLcLtkrFlags;
extern G__linked_taginfo G__ReconCintLN_ReconHeadercLcLcalFlags;
extern G__linked_taginfo G__ReconCintLN_Recon;
extern G__linked_taginfo G__ReconCintLN_TkrLocator;
extern G__linked_taginfo G__ReconCintLN_TkrHit;
extern G__linked_taginfo G__ReconCintLN_TkrLocatorcLcLLocationType;
extern G__linked_taginfo G__ReconCintLN_TkrSiClustercLcLTKRAxes;
extern G__linked_taginfo G__ReconCintLN_TkrTrackcLcLParticleType;

/* STUB derived class for protected member access */
