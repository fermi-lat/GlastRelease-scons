
#ifndef GlastMcTask_h
#define GlastMcTask_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TFolder.h>
#include <TTask.h>
#include <TCut.h>
#if !defined(__CINT__)
#include <iostream>
#else  // for interactive use
#include "Riostream.h"
#endif

class RequestedMcData : public TObject {
public:
    RequestedMcData() { m_dataWord = 0; };
    virtual ~RequestedMcData() {};

    typedef enum { 
        INTHITCOL = 1,
        POSHITCOL = 2,
        PARTICLECOL = 4,
        ALL = 7
    } McData;

   // Unpack the requested dataWord to determine what types are required
   bool McIntHitData() const { return m_dataWord & INTHITCOL; };
   bool McPosHitData() const { return m_dataWord & POSHITCOL; };
   bool McPartData() const { return m_dataWord & PARTICLECOL; };

   void RegisterData(RequestedMcData::McData d){ m_dataWord |= d;};

private:
    unsigned int m_dataWord;

    ClassDef(RequestedMcData,0)
};

class GlastMcTask : public TTask {
public:

   GlastMcTask():TTask() { };
   GlastMcTask(const char *name,const char *title);
   virtual ~GlastMcTask() { ; }
   void Exec(Option_t *option="");
   void CreateHistograms();
   virtual void CreateHistograms(TFolder* folder){};

   void RegisterData(RequestedMcData::McData d);

protected:
   TFolder *glastFolder, *histoFolder, *taskFolder;

private:

 ClassDef(GlastMcTask,0)
};


#endif

