
#ifndef AcdReconTask_h
#define AcdReconTask_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTask.h>
#include <TCut.h>
#include "GlastReconTask.h"
#if !defined(__CINT__)
#include <iostream>
#else  // for interactive use
#include "Riostream.h"
#endif


class AcdReconTask : public GlastReconTask {
public:
   AcdReconTask();
   AcdReconTask(const char *name,const char *title);
   virtual ~AcdReconTask() { ; }
   void Exec(Option_t *option="");
   void CreateHistograms(TFolder *folder);

private:

 ClassDef(AcdReconTask,0)
};

#endif

