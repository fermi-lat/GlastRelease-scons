
#ifndef AcdDigiTask_cxx
#define AcdDigiTask_cxx

#include "AcdDigiTask.h"
#include "TH1.h"
#include "digiRootData/DigiEvent.h"
#include "TFolder.h"

ClassImp(AcdDigiTask)

AcdDigiTask::AcdDigiTask() :GlastDigiTask() {
    GlastDigiTask::RegisterData(RequestedDigiData::ACD);
}

AcdDigiTask::AcdDigiTask(const char *name,const char *title):GlastDigiTask(name,title)
{
    GlastDigiTask::RegisterData(RequestedDigiData::ACD);
}


void AcdDigiTask::CreateHistograms(TFolder *folder) {

    TH1F *PMT0001 = new TH1F("PMT0001", "PMT0001", 1000, 0, 1000);
    PMT0001->SetDirectory(0);
    folder->Add(PMT0001);

}

void AcdDigiTask::Exec(Option_t *option) {
    TH1F* pmt0001 = (TH1F*)gROOT->FindObjectAny("PMT0001");
    if (!pmt0001) {
        std::cout << "PMT0001 Histogram Not Found!" << std::endl;
        return;
    }
    TFolder *digiFolder = (TFolder*)gROOT->FindObjectAny("/glast/digi"); 
    TClonesArray *m_acdDigiCol = (TClonesArray*)digiFolder->FindObject("AcdDigis");
    if (!m_acdDigiCol) {
        std::cout << "Could not retrieve the digiEvent" << std::endl;
        return;
    }
    
    AcdId id0001(0,0,0,1);
    TIter acdDigiIter(m_acdDigiCol);
    AcdDigi *acdDigiItem = 0;
    while (acdDigiItem = (AcdDigi*)acdDigiIter.Next()) {
        AcdId id = acdDigiItem->getId();
        if (id.getId() == id0001.getId()) {
            pmt0001->Fill(acdDigiItem->getPulseHeight(AcdDigi::A));
        }
    }
}

#endif

