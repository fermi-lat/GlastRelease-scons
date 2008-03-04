/** 
*@class takedata_mt
*
*@brief Root Tree Ananlysis Trigger time structure of a run
*
* This class is for Trigger Root Tree Ananlysis of the run time structure
*
* Version 0.1 30-Jun-2006 Martin Kocian Initial version
*/

#ifndef takedata_mt_cxx
#define takedata_mt_cxx 1

#include "calibGenTRG/treqACD.h"
#include "TestReport.h"
#include "eval.h"
#include "tackcuts.h"

#include <TF1.h>
#include <TCanvas.h>
#include "calibGenTRG/RootTreeAnalysis.h"
#include "TProfile.h"

#include "CalGeom.h"

#include <list>
#include <time.h>

using namespace CLHEP;
using namespace CalUtil;

ClassImp(treqACD)
  
// 
// Initialization
//
treqACD::treqACD() 
{
}

void treqACD::setParameters(int tkrdelay, int acddelay){
  m_tkrdelay=tkrdelay;
  m_acddelay=acddelay;
  std::cout<<"TKR TREQ delay: "<<m_tkrdelay<<std::endl;
  std::cout<<"ACD TREQ delay: "<<m_tkrdelay<<std::endl;
}

void treqACD::inithistos(){

  //  char histname[128];
    //char histtitle[128];
    //sprintf(histname,"ACD_adchist");
    //sprintf(histtitle,"ACD ADC histogram");
    //m_acdadchist=new TH1D(histname,histtitle,100,0,3);
    //m_acdadchist->GetXaxis()->SetTitle("MIPS");
  m_allevents=new TH1D("allevents","All Events",32,-.5,31.5);
}  

void treqACD::savehistos(char* filename){
  gROOT->SetStyle("Plain");
  TCanvas c1;
  m_allevents->Draw();
  c1.SaveAs("allevents.gif");
  strcpy(m_filename,filename);
  TFile f(filename,"recreate");
  //  TF1* langau=new TF1(LangauFun::getLangauDAC());
    //langau->SetParameters(.1,1,30000,.15,1000);
    //m_acdadchist->Fit(langau,"em","",0.1,3);
    //m_mp=langau->GetParameter(1);
    //m_mp_err=langau->GetParError(1);
    //m_lanwid=langau->GetParameter(0);
    //m_lanwid_err=langau->GetParError(0);
    //m_gaussig=langau->GetParameter(3);
    //m_gaussig_err=langau->GetParError(3);
    //m_chi2=langau->GetChisquare();
  m_allevents->Write();
  f.Close();
}
  
void treqACD::writeoutresults(char* reportname){
  char name[128];
  TestReport r(reportname);
  r.newheadline("Test identification");
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  sprintf(name,"%s",asctime(timeinfo));
  r.additem("Date", name);
  r.additem("Application", "vetoeff");
  r.additem("Description", "Analysis of a veto efficiency run");
  r.additem("Type of analysis","Offline");
  sprintf (name,"%d",m_gid);
  r.additem("Run ground id",name);
  sprintf (name,"%d",m_nev);
  r.additem("Number of events",name);
  if (m_status==0)r.addstatus("Passed");
  else r.addstatus("Failed");
  r.newheadline("Results");
  r.addimage("allevents.gif","ACD Conditions arrival time for all events.");
  r.writereport();
}
//
// Main event loop
//
void treqACD::Go(Long64_t numEvents)
{        
    // For this dump example, activate read in of all tree branches

    ActivateAllBranches(); 

    // To speed up read in, can choose to disable some branches 
    digiTree->SetBranchStatus("*",1);  // enable all branches
    //reconTree->SetBranchStatus("*",1);  // disable recon branches
    
       
    // determine how many events to process
    Long64_t nentries = GetEntries();
    std::cout << "\nNum Events in File is: " << nentries << std::endl;

    Int_t nEvtMax = TMath::Min(numEvents+m_StartEvent,nentries);
    if (m_StartEvent == nentries) {
        std::cout << " all events in file read" << std::endl;
        return;
    }
    if (nentries <= 0) return;
    if (mc) mc->Clear();
    if (evt) evt->Clear();
    if (rec) rec->Clear();
    
    Int_t nb = GetEvent(m_StartEvent);
    m_gid=evt->getMetaEvent().run().id();
    m_gid+=77000000;
    std::cout<<"gid "<<m_gid<<std::endl;
    m_starttime=evt->getTimeStamp( );
    m_nev=nEvtMax-m_StartEvent;

    //}
    // BEGINNING OF EVENT LOOP
    for (Long64_t ievent=m_StartEvent; ievent<nEvtMax; ievent++) {
        
        if (mc) mc->Clear();
        if (evt) evt->Clear();
        if (rec) rec->Clear();
        
        Int_t nb = GetEvent(ievent);

        if(ievent%1000==0) 
            std::cout << "** Processing Event " << ievent << std::endl;  
	const Gem gem=evt->getGem();
	if(gem.getConditionSummary()==3){  // TKR && ROI
	  UShort_t condarracd=gem.getCondArrTime().roi();
	  m_allevents->Fill((double)condarracd);
	}
	/*
	if (rec){

	  AcdRecon *acdRec = rec->getAcdRecon();
	  if (acdRec){
	    int nhit=acdRec->nAcdHit();
	    for (int i=0;i<nhit;i++){
	      const AcdHit* acdHit=acdRec->getAcdHit(i);
	      bool intersect=false;
	      double pathlength=10;
	      for (unsigned j=0;j<acdRec->nAcdIntersections();j++){
		if (acdHit->getId()==acdRec->getAcdTkrIntersection(j)->getTileId()){
		  pathlength=acdRec->getAcdTkrIntersection(j)->getPathLengthInTile()/10.;
		  intersect=true;
		  break;
		}
	      }
	      if(intersect){
		//std::cout<<"found intersection"<<std::endl;
		double acdenergy=(double)acdHit->getMips()/pathlength;
		
		m_acdadchist->Fill(acdenergy);
		//		if(acdenergy==0)std::cout<<"0 energy in tile "<<acdHit->getId().getId()<<std::endl;
	      }
	    }	
	  }  
	}
	*/
    }  // end analysis code in event loop
    m_status=1;
    

std::cout << "**** End of event loop **** " << std::endl;


}



#endif
