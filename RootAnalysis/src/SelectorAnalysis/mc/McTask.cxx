
#ifndef McTask_cxx
#define McTask_cxx

#include "McTask.h"
#include "TH1.h"
#include "mcRootData/McEvent.h"
#include "TFolder.h"

ClassImp(McTask)

McTask::McTask() :GlastMcTask() {
    GlastMcTask::RegisterData(RequestedMcData::POSHITCOL);
    GlastMcTask::RegisterData(RequestedMcData::INTHITCOL);
    GlastMcTask::RegisterData(RequestedMcData::PARTICLECOL);
}

McTask::McTask(const char *name,const char *title):GlastMcTask(name,title)
{
    GlastMcTask::RegisterData(RequestedMcData::POSHITCOL);
    GlastMcTask::RegisterData(RequestedMcData::INTHITCOL);
    GlastMcTask::RegisterData(RequestedMcData::PARTICLECOL);
}


void McTask::CreateHistograms(TFolder *folder) {

    TH1F *McPartId = new TH1F("McPartId", "McPartId", 20, 0.5, 20.5);
    McPartId->SetDirectory(0);
    folder->Add(McPartId);

}

void McTask::Exec(Option_t *option) {
    TH1F* McPartId = (TH1F*)gROOT->FindObjectAny("McPartId");
    if (!McPartId) {
        std::cout << "McPartId Histogram Not Found!" << std::endl;
        return;
    }
    TFolder *mcFolder = (TFolder*)gROOT->FindObjectAny("/glast/mc"); 
    McEvent *m_mcEvent = (McEvent*)mcFolder->FindObject("McEvent");
    if (!m_mcEvent) {
      std::cout << "could not retrieve mcEvent" << std::endl;
      return;
    }
//    TObjArray *m_mcPartCol = (TObjArray*)mcFolder->FindObject("McParticle");
  //  if (!m_mcPartCol) {
  //      std::cout << "Could not retrieve the mcPartCol" << std::endl;
  //      return;
   // }
    
    TObjArray *m_mcPartCol = m_mcEvent->getMcParticleCol();
    if (!m_mcPartCol) return;
    TIter mcPartIter(m_mcPartCol);
    McParticle *mcPart = 0;
    while (mcPart = (McParticle*)mcPartIter.Next()) {
        Float_t id = mcPart->getParticleId();
        McPartId->Fill(id);
    }
}

#endif

