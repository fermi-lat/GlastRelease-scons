
#ifndef GlastDigiTask_cxx
#define GlastDigiTask_cxx

#include "GlastDigiTask.h"
#include "TH1.h"
#include "TFolder.h"

ClassImp(RequestedDigiData)
ClassImp(GlastDigiTask)

GlastDigiTask::GlastDigiTask(const char *name,const char *title) : TTask(name, title)
{
    glastFolder = (TFolder*)gROOT->FindObjectAny("/glast");
    if (!glastFolder)
        glastFolder = gROOT->GetRootFolder()->AddFolder("glast","glast top level folders");
    taskFolder = (TFolder*)glastFolder->FindObjectAny("/glast/tasks");
    if (!taskFolder)
        taskFolder = glastFolder->AddFolder("tasks", "tasks folder");
}

void GlastDigiTask::RegisterData(RequestedDigiData::DigiData d) {
    RequestedDigiData *requested = (RequestedDigiData*)taskFolder->FindObjectAny("RequestedDigiData");
    if (!requested) {
        RequestedDigiData *DigiData = new RequestedDigiData();
        DigiData->RegisterData(d);
        taskFolder->Add(DigiData);
    } else {
        requested->RegisterData(d);
    }
    return;
}

void GlastDigiTask::CreateHistograms() {
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
      ((GlastDigiTask*)task)->CreateHistograms(histoFolder);
    }     
}

void GlastDigiTask::Exec(Option_t *option) {


}




#endif

