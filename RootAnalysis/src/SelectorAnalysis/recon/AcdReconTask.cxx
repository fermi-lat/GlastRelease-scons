
#ifndef AcdReconTask_cxx
#define AcdReconTask_cxx

#include "AcdReconTask.h"
#include "TH1.h"
#include "reconRootData/ReconEvent.h"
#include "TFolder.h"

ClassImp(AcdReconTask)

AcdReconTask::AcdReconTask() :GlastReconTask() {
    GlastReconTask::RegisterData(RequestedReconData::ACD);
}

AcdReconTask::AcdReconTask(const char *name,const char *title):GlastReconTask(name,title)
{
    GlastReconTask::RegisterData(RequestedReconData::ACD);
}


void AcdReconTask::CreateHistograms(TFolder *folder) {

    TH1F *AcdActDist = new TH1F("AcdActDist", "AcdActDist", 1000, 0, 1000);
    AcdActDist->SetDirectory(0);
    folder->Add(AcdActDist);

}

void AcdReconTask::Exec(Option_t *option) {
    TH1F* AcdActDist = (TH1F*)gROOT->FindObjectAny("AcdActDist");
    if (!AcdActDist) {
        std::cout << "AcdActDist Histogram Not Found!" << std::endl;
        return;
    }
    TFolder *reconFolder = (TFolder*)gROOT->FindObjectAny("/glast/recon"); 
    ReconEvent *rec = (ReconEvent*)reconFolder->FindObject("ReconEvent");
    if (!rec) return;
    AcdRecon *m_acdRecon = rec->getAcdRecon();
    //AcdRecon *m_acdRecon = (AcdRecon*)reconFolder->FindObject("AcdRecon");
    if (!m_acdRecon) {
        std::cout << "Could not retrieve the AcdRecon" << std::endl;
        return;
    }
    
    AcdActDist->Fill(m_acdRecon->getActiveDist());
}

#endif

