
#ifndef GlastMcTask_cxx
#define GlastMcTask_cxx

#include "GlastMcTask.h"
#include "TH1.h"
#include "TFolder.h"

ClassImp(RequestedMcData)
ClassImp(GlastMcTask)

GlastMcTask::GlastMcTask(const char *name,const char *title) : TTask(name, title)
{
    glastFolder = (TFolder*)gROOT->FindObjectAny("/glast");
    if (!glastFolder)
        glastFolder = gROOT->GetRootFolder()->AddFolder("glast","glast top level folders");
    taskFolder = (TFolder*)glastFolder->FindObjectAny("/glast/tasks");
    if (!taskFolder)
        taskFolder = glastFolder->AddFolder("tasks", "tasks folder");
}

void GlastMcTask::RegisterData(RequestedMcData::McData d) {
    RequestedMcData *requested = (RequestedMcData*)taskFolder->FindObjectAny("RequestedMcData");
    if (!requested) {
        RequestedMcData *McData = new RequestedMcData();
        McData->RegisterData(d);
        taskFolder->Add(McData);
    } else {
        requested->RegisterData(d);
    }
    return;
}

void GlastMcTask::CreateHistograms() {
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
      ((GlastMcTask*)task)->CreateHistograms(histoFolder);
    }     
}

void GlastMcTask::Exec(Option_t *option) {


}




#endif

