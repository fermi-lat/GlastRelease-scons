
#ifndef GlastReconTask_cxx
#define GlastReconTask_cxx

#include "GlastReconTask.h"
#include "TH1.h"
#include "TFolder.h"

ClassImp(RequestedReconData)
ClassImp(GlastReconTask)

GlastReconTask::GlastReconTask(const char *name,const char *title) : TTask(name, title)
{
    glastFolder = (TFolder*)gROOT->FindObjectAny("/glast");
    if (!glastFolder)
        glastFolder = gROOT->GetRootFolder()->AddFolder("glast","glast top level folders");
    taskFolder = (TFolder*)glastFolder->FindObjectAny("/glast/tasks");
    if (!taskFolder)
        taskFolder = glastFolder->AddFolder("tasks", "tasks folder");
}

void GlastReconTask::RegisterData(RequestedReconData::ReconData d) {
    RequestedReconData *requested = (RequestedReconData*)taskFolder->FindObjectAny("RequestedReconData");
    if (!requested) {
        RequestedReconData *ReconData = new RequestedReconData();
        ReconData->RegisterData(d);
        taskFolder->Add(ReconData);
    } else {
        requested->RegisterData(d);
    }
    return;
}

void GlastReconTask::CreateHistograms() {
    histoFolder = (TFolder*)gROOT->FindObjectAny("/glast/histograms");
    if (!histoFolder) {
        glastFolder = (TFolder*)gROOT->FindObjectAny("/glast");
        if (!glastFolder) 
            glastFolder = gROOT->GetRootFolder()->AddFolder("glast","glast top level folders");
        glastFolder->AddFolder("histograms", "histogram folder");
    }
    TIter next(fTasks);
    TTask *task;
    while((task=(TTask*)next())) {
      if (!task->IsActive()) continue;
      ((GlastReconTask*)task)->CreateHistograms(histoFolder);
    }     
}

void GlastReconTask::Exec(Option_t *option) {


}




#endif

