/** 
*@class afterglow
*
*@brief Root Tree Ananlysis 
*
* This class is for Trigger Root Tree Ananlysis of the run time structure
*
* Version 0.1 30-Jun-2006 Martin Kocian Initial version
*/

#ifndef afterglow_cxx
#define afterglow_cxx 1

#include "calibGenTRG/afterglow.h"
#include "eval.h"
#include "tackcuts.h"

#include "calibGenTRG/RootTreeAnalysis.h"
#include "TProfile.h"
#include "CalGeom.h"

#include <list>
#include <set>

using namespace CLHEP;
using namespace CalUtil;

ClassImp(afterglow)
  
// 
// Initialization
//
afterglow::afterglow() 
{
}

void afterglow::inithistos(){

  m_deltafrac=new TH2D("deltafrac","Fraction of TKR hits left",500,0,5000,100,0,1);
  m_deltafrac->SetMarkerStyle(22);
  m_deltalayerfrac=new TH2D("deltalayerfrac","Fraction of TKR layers left",500,0,5000,100,0,1);
  m_deltalayerfrac->SetMarkerStyle(22);
  m_deltalayerfraceng=new TH2D("deltalayerfraceng","Fraction of TKR layers left",500,0,5000,100,0,1);
  m_deltalayerfracengnorm=new TH2D("deltalayerfracengnorm","Fraction of TKR layers left",500,0,5000,100,0,1);
  m_hits=new TH2D("hits","Number of hits in previous event",500,0,5000,200,0,1000);
  m_hits->SetMarkerStyle(22);
  m_layers=new TH2D("layers","Number of layers in previous event",500,0,5000,40,0,40);
  m_layers->SetMarkerStyle(22);
  m_delta=new TH1D("delta","delta event time for periodic triggers",500,0,5000);

}  

void afterglow::savehistos(char* filename){
  strcpy(m_filename,filename);
  TFile f(filename,"recreate");
  m_deltafrac->Write();
  m_delta->Write();
  m_deltalayerfrac->Write();
  m_deltalayerfraceng->Divide(m_deltalayerfracengnorm);
  m_deltalayerfraceng->Write();
  m_hits->Write();
  m_layers->Write();
  f.Close();
}
  
//
// Main event loop
//
void afterglow::Go(Long64_t numEvents)
{        
    // For this dump example, activate read in of all tree branches

    //ActivateAllBranches(); 

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
    
    //Int_t nb = GetEvent(m_StartEvent);
    //m_gid=evt->getMetaEvent().run().id();
    //std::cout<<"gid "<<m_gid<<std::endl;
    m_nev=nEvtMax-m_StartEvent;
    bool cno=false;
    std::vector<unsigned> oldhits;
    //std::set<unsigned> oldlayers;
    double oldeng=0;
    //}
    // BEGINNING OF EVENT LOOP
    for (Long64_t ievent=m_StartEvent; ievent<nEvtMax; ievent++) {
        
        if (mc) mc->Clear();
        if (evt) evt->Clear();
        if (rec) rec->Clear();
        
        GetEvent(ievent);

        if(ievent%10000==0) 
            std::cout << "** Processing Event " << ievent << std::endl;  
	UShort_t cond=evt->getGem().getConditionSummary();
	UShort_t delta=evt->getGem().getDeltaEventTime();
	if ((cond&32)|| (cond&3)==3){
	  //if (cond&48){
	  std::vector<unsigned> hits;
	  //std::set<unsigned> layers;
	  const TObjArray* tkrDigiCol = evt->getTkrDigiCol();
	  if (tkrDigiCol){
	    TIter tkrIter(tkrDigiCol);
	    TkrDigi *tkr = 0;
	    while ((tkr = (TkrDigi*)tkrIter.Next())) {
	      Int_t tower   = tkr->getTower().id();
	      Int_t bilayer = tkr->getBilayer();
	      Int_t view = tkr->getView();
	      //layers.insert(view*10000+bilayer*100000+tower*10000000);
	      for (unsigned int i=0;i<tkr->getNumHits( );i++){
		unsigned strip=tkr->getStrip(i);
		hits.push_back(strip+view*10000+bilayer*100000+tower*10000000);
	      }
	    }
	  }  
	  double eng=0;
	  if(rec){
	    TObjArray *xtalRecCol = rec->getCalRecon()->getCalXtalRecCol();
	    TIter xtalIter(xtalRecCol);
	    CalXtalRecData *xtal = 0;
	    std::set<unsigned> callayers;
	    if (xtalRecCol){
	      if(xtalRecCol->GetEntries()>0 && xtalRecCol->GetEntries()<10){
		while ((xtal = (CalXtalRecData*)xtalIter.Next())) {
		  eng+=xtal->getEnergy();
		  callayers.insert(xtal->getPackedId().getLayer());
		}
		eng/=callayers.size();
	      }
	    }
	  }
	  if (!(cond&32)){
	    if(hits.size()<80)cno=true;
	    else cno=false;
	    oldhits=hits;
	    //oldlayers=layers;
	    oldeng=eng;
	  }else if (cond==32){
	    if (delta<5000)m_delta->Fill(delta);
	    if (cno==true){
	      if (oldhits.size()>30){
		unsigned int found=0;
		//unsigned int foundlayers=0;
		std::set<unsigned> layers;
		std::set<unsigned> oldlayers;
		for (unsigned int i=0;i<oldhits.size();i++){
		  oldlayers.insert(oldhits[i]/10000);
		  if (std::find(hits.begin(),hits.end(),oldhits[i])!=hits.end()){
		    found++;
		    layers.insert(oldhits[i]/10000);
		  }
		}
		//std::set<unsigned>::iterator it;
		//for (it=oldlayers.begin();it!=oldlayers.end();it++){
		//  if (std::find(layers.begin(),layers.end(),(*it))!=layers.end())foundlayers++;
		//}
		double frac=(double)found/double(oldhits.size());
		double layerfrac=(double)layers.size()/double(oldlayers.size());
		if (delta<5000){
		  m_deltafrac->Fill(delta,frac);
		  m_hits->Fill(delta,oldhits.size());
		  m_deltalayerfrac->Fill(delta,layerfrac);
		  if (delta<1500&&layerfrac<0.7)std::cout<<"Event "<<ievent<<" delta "<<delta<<" frac "<<layerfrac<<std::endl;
		  if (delta>2200&&layerfrac>0.6)std::cout<<"Event "<<ievent<<" delta "<<delta<<" frac "<<layerfrac<<std::endl;
		  if (oldeng>0&&oldeng<1000){
		    m_deltalayerfraceng->Fill(delta,layerfrac,oldeng);
		    m_deltalayerfracengnorm->Fill(delta,layerfrac);
		  }
		  m_layers->Fill(delta,oldlayers.size());
		}
	      }
	      cno=false;
	    }
	  }else{
	    cno=false;
	  }
	}else{
	  cno=false;
	}
	    
	      
    }  // end analysis code in event loop
    

std::cout << "**** End of event loop **** " << std::endl;

// GEM time counter summary
//igemUtil->dumpTimeCounter();

// Endjob histograms  
//AcdHend();

}



#endif
