
#ifndef GlastReconTask_h
#define GlastReconTask_h

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

class RequestedReconData : public TObject {
public:
    RequestedReconData() { m_dataWord = 0; };
    virtual ~RequestedReconData() {};

    typedef enum { 
        ACD = 1,
        CAL = 2,
        TKR = 4,
        ALL = 7
    } ReconData;

   // Unpack the requested dataWord to determine what types are required
   bool AcdData() const { return m_dataWord & ACD; };
   bool CalData() const { return m_dataWord & CAL; };
   bool TkrData() const { return m_dataWord & TKR; };

   void RegisterData(RequestedReconData::ReconData d){ m_dataWord |= d;};

private:
    unsigned int m_dataWord;

    ClassDef(RequestedReconData,0)
};

class GlastReconTask : public TTask {
public:

   GlastReconTask():TTask() { };
   GlastReconTask(const char *name,const char *title);
   virtual ~GlastReconTask() { ; }
   void Exec(Option_t *option="");
   void CreateHistograms();
   virtual void CreateHistograms(TFolder* folder){};

   void RegisterData(RequestedReconData::ReconData d);

protected:
   TFolder *glastFolder, *histoFolder, *taskFolder;

private:

 ClassDef(GlastReconTask,0)
};


#endif

