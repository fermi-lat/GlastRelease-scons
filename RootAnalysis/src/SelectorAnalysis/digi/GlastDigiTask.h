
#ifndef GlastDigiTask_h
#define GlastDigiTask_h

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

class RequestedDigiData : public TObject {
public:
    RequestedDigiData() { m_dataWord = 0; };
    virtual ~RequestedDigiData() {};

    typedef enum { 
        ACD = 1,
        CAL = 2,
        TKR = 4,
        EBF = 8,
        HEADER = 16,
        L1T = 32,
        ALL = 53
    } DigiData;

   // Unpack the requested dataWord to determine what types are required
   bool AcdData() const { return m_dataWord & ACD; };
   bool CalData() const { return m_dataWord & CAL; };
   bool TkrData() const { return m_dataWord & TKR; };
   bool EbfData() const { return m_dataWord & EBF; };
   bool HeaderData() const { return m_dataWord & HEADER; };
   bool L1TData() const { return m_dataWord & L1T; };

   void RegisterData(RequestedDigiData::DigiData d){ m_dataWord |= d;};

private:
    unsigned int m_dataWord;

    ClassDef(RequestedDigiData,0)
};

class GlastDigiTask : public TTask {
public:

   GlastDigiTask():TTask() { };
   GlastDigiTask(const char *name,const char *title);
   virtual ~GlastDigiTask() { ; }
   void Exec(Option_t *option="");
   void CreateHistograms();
   virtual void CreateHistograms(TFolder* folder){};

   void RegisterData(RequestedDigiData::DigiData d);

protected:
   TFolder *glastFolder, *histoFolder, *taskFolder;

private:

 ClassDef(GlastDigiTask,0)
};


#endif

