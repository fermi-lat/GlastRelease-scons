
#ifndef McTask_h
#define McTask_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTask.h>
#include <TCut.h>
#include "GlastMcTask.h"
#if !defined(__CINT__)
#include <iostream>
#else  // for interactive use
#include "Riostream.h"
#endif


class McTask : public GlastMcTask {
public:
   McTask();
   McTask(const char *name,const char *title);
   virtual ~McTask() { ; }
   void Exec(Option_t *option="");
   void CreateHistograms(TFolder *folder);

private:

 ClassDef(McTask,0)
};

#endif

