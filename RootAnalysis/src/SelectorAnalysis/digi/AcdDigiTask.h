
#ifndef AcdDigiTask_h
#define AcdDigiTask_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTask.h>
#include <TCut.h>
#include "GlastDigiTask.h"
#if !defined(__CINT__)
#include <iostream>
#else  // for interactive use
#include "Riostream.h"
#endif


class AcdDigiTask : public GlastDigiTask {
public:
   AcdDigiTask();
   AcdDigiTask(const char *name,const char *title);
   virtual ~AcdDigiTask() { ; }
   void Exec(Option_t *option="");
   void CreateHistograms(TFolder *folder);

private:

 ClassDef(AcdDigiTask,0)
};

#endif

