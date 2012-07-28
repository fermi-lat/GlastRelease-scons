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

#include "calibGenTRG/mippeak.h"
#include "eval.h"
#include "tackcuts.h"

#include <TF1.h>
#include "calibGenTRG/RootTreeAnalysis.h"
#include "LangauFun.h"
#include "TProfile.h"
#include "CalGeom.h"

#include <list>

using namespace CLHEP;
using namespace CalUtil;

ClassImp(mippeak)
  
// 
// Initialization
//
mippeak::mippeak() 
{
}

void mippeak::inithistos(){

  char histname[128];
  char histtitle[128];
  sprintf(histname,"ACD_adchist");
  sprintf(histtitle,"ACD ADC histogram");
  m_acdadchist=new TH1D(histname,histtitle,100,0,3);
  m_acdadchist->GetXaxis()->SetTitle("MIPS");
}  

void mippeak::savehistos(char* filename){
  strcpy(m_filename,filename);
  TFile f(filename,"recreate");
  TF1* langau=new TF1(LangauFun::getLangauDAC());
  langau->SetParameters(.1,1,30000,.15,1000);
  m_acdadchist->Fit(langau,"em","",0.1,3);
  m_mp=langau->GetParameter(1);
  m_mp_err=langau->GetParError(1);
  m_lanwid=langau->GetParameter(0);
  m_lanwid_err=langau->GetParError(0);
  m_gaussig=langau->GetParameter(3);
  m_gaussig_err=langau->GetParError(3);
  m_chi2=langau->GetChisquare();
  m_acdadchist->Write();
  f.Close();
}
  
void mippeak::writeoutresults(char* rf){
  ofstream res(rf);
  res<<m_gid<<std::endl;
  res<<m_starttime/3600./24.<<std::endl;
  res<<m_chi2<<std::endl;
  res<<m_mp<<std::endl;
  res<<m_mp_err<<std::endl;
  res<<m_lanwid<<std::endl;
  res<<m_lanwid_err<<std::endl;
  res<<m_gaussig<<std::endl;
  res<<m_gaussig_err<<std::endl;

  
}



//
// Main event loop
//
void mippeak::Go(Long64_t numEvents)
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
    
    GetEvent(m_StartEvent);
    m_gid=evt->getMetaEvent().run().id();
    std::cout<<"gid "<<m_gid<<std::endl;
    m_starttime=evt->getTimeStamp( );
    m_nev=nEvtMax-m_StartEvent;

    //}
    // BEGINNING OF EVENT LOOP
    for (Long64_t ievent=m_StartEvent; ievent<nEvtMax; ievent++) {
        
        if (mc) mc->Clear();
        if (evt) evt->Clear();
        if (rec) rec->Clear();
        
        GetEvent(ievent);

        if(ievent%1000==0) 
            std::cout << "** Processing Event " << ievent << std::endl;  
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
    }  // end analysis code in event loop
    

std::cout << "**** End of event loop **** " << std::endl;

// GEM time counter summary
//igemUtil->dumpTimeCounter();

// Endjob histograms  
//AcdHend();

}



#endif
