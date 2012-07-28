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
  std::cout<<"ACD TREQ delay: "<<m_acddelay<<std::endl;
}

void treqACD::inithistos(){

  char histname[128];
  char histtitle[128];
    //sprintf(histname,"ACD_adchist");
    //sprintf(histtitle,"ACD ADC histogram");
    //m_acdadchist=new TH1D(histname,histtitle,100,0,3);
    //m_acdadchist->GetXaxis()->SetTitle("MIPS");
  m_allevents=new TH1D("allevents","All Events",32,-.5,31.5);
  m_allevents->GetXaxis()->SetTitle("ACD conditions arrival time (ticks)");
  m_intersect=new TH1D("intersect","Events where Track intersects tile",32,-.5,31.5);
  m_intersect->GetXaxis()->SetTitle("ACD conditions arrival time (ticks)");
  m_singletile=new TH1D("singletile","Single Tile Events",32,-.5,31.5);
  m_singletile->GetXaxis()->SetTitle("ACD conditions arrival time (ticks)");
  m_singletower=new TH1D("singletower","Single Tower Events",32,-.5,31.5);
  m_singletower->GetXaxis()->SetTitle("ACD conditions arrival time (ticks)");
  m_singletowertile=new TH1D("singletowertile","Single Tower and Tile Events",32,-.5,31.5);
  m_singletowertile->GetXaxis()->SetTitle("ACD conditions arrival time (ticks)");
  m_cno=new TH1D("cno","CNO",32,-.5,31.5);
  m_cno->GetXaxis()->SetTitle("ACD conditions arrival time (ticks)");
  for (int i=0;i<16;i++){
    sprintf(histname,"tower_%d",i);
    sprintf(histtitle,"Tower %d",i);
    m_bytower[i]=new TH1D(histname,histtitle,32,-.5,31.5);
    m_bytower[i]->GetXaxis()->SetTitle("ACD conditions arrival time (ticks)");
  }
  for (int i=0;i<108;i++){
    unsigned int tile=AcdId::tileFromIndex(i);
    sprintf(histname,"tile_%03d",tile);
    sprintf(histtitle,"Tile %03d",tile);
    m_bytile[i]=new TH1D(histname,histtitle,32,-.5,31.5);
    m_bytile[i]->GetXaxis()->SetTitle("ACD conditions arrival time (ticks)");
  }
}  


void treqACD::fit(TH1D* histo,double& mean, double &emean, double& sigma, double& esigma, int& nevents){
  histo->Fit("gaus");
  mean=histo->GetFunction("gaus")->GetParameter(1);
  emean=histo->GetFunction("gaus")->GetParError(1);
  sigma=histo->GetFunction("gaus")->GetParameter(2);
  esigma=histo->GetFunction("gaus")->GetParError(2);
  nevents=(int)histo->GetEntries();
}

void treqACD::writeoutresults(const char* reportname, const char* filename){
  m_status=0;
  char name[128];
  char title[128];
  char* line[7];
  for (int j=0;j<7;j++)line[j]=new char[64];
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(10);
  gStyle->SetOptFit(11);
  TCanvas c1;
  TFile f(filename,"recreate");
  TestReport r(const_cast<char*>(reportname));
  r.newheadline("Test identification");
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  sprintf(name,"%s",asctime(timeinfo));
  r.additem("Date", name);
  r.additem("Application", "treqACD");
  r.additem("Description", "Analysis of a Treq_TkrAcd run");
  r.additem("Type of analysis","Offline");
  sprintf (name,"%d",m_gid);
  r.additem("Run ground id",name);
  sprintf (name,"%d",m_nev);
  r.additem("Number of events",name);
  sprintf (name,"%d",m_tkrdelay);
  r.additem("TKR delay",name);
  sprintf (name,"%d",m_acddelay);
  r.additem("ACD delay",name);

  r.newheadline("General Results");
  r.addimage("allevents.gif");
  r.addimage("cno.gif");
  char* restable[]={"Cut","Number of events","Mean +- error","Sigma +- error","Optimal delay difference"};
  double mean,emean,sigma,esigma;
  int nentries, bestdelay, delay;
  r.starttable(restable,5);
  fit(m_allevents,mean,emean,sigma,esigma,nentries);
  m_mean=mean;
  m_sigma=sigma;
  bestdelay=(int)(mean+.5);
  delay=m_acddelay-m_tkrdelay-bestdelay;
  r.linktext(line[0],"All Events", "allevents.gif");
  sprintf(line[1],"%d",nentries);
  if(nentries<1000)markError(line[1],&r);
  sprintf(line[2],"%.2f +- %.2f",mean,emean);
  sprintf(line[3],"%.2f +- %.2f",sigma,esigma);
  sprintf(line[4],"%d",delay);
  r.addtableline(line,5);
  c1.SaveAs("allevents.gif");
  m_allevents->Write();
  
  fit(m_intersect,mean,emean,sigma,esigma,nentries);
  bestdelay=(int)(mean+.5);
  delay=m_acddelay-m_tkrdelay-bestdelay;
  r.linktext(line[0],"Events with TKR ACD Intersection", "intersect.gif");
  sprintf(line[1],"%d",nentries);
  if(nentries<1000)markError(line[1],&r);
  sprintf(line[2],"%.2f +- %.2f",mean,emean);
  if(fabs(mean-m_mean)>1)markError(line[2],&r);
  sprintf(line[3],"%.2f +- %.2f",sigma,esigma);
  if(fabs(sigma-m_sigma)>2)markError(line[3],&r);
  sprintf(line[4],"%d",delay);
  r.addtableline(line,5);
  c1.SaveAs("intersect.gif");
  m_intersect->Write();

  fit(m_singletile,mean,emean,sigma,esigma,nentries);
  bestdelay=(int)(mean+.5);
  delay=m_acddelay-m_tkrdelay-bestdelay;
  r.linktext(line[0],"Single Tile Events", "singletile.gif");
  sprintf(line[1],"%d",nentries);
  if(nentries<1000)markError(line[1],&r);
  sprintf(line[2],"%.2f +- %.2f",mean,emean);
  if(fabs(mean-m_mean)>1)markError(line[2],&r);
  sprintf(line[3],"%.2f +- %.2f",sigma,esigma);
  if(fabs(sigma-m_sigma)>2)markError(line[3],&r);
  sprintf(line[4],"%d",delay);
  r.addtableline(line,5);
  c1.SaveAs("singletile.gif");
  m_singletile->Write();
  
  fit(m_singletower,mean,emean,sigma,esigma,nentries);
  bestdelay=(int)(mean+.5);
  delay=m_acddelay-m_tkrdelay-bestdelay;
  r.linktext(line[0],"Single Tower Events", "singletower.gif");
  sprintf(line[1],"%d",nentries);
  if(nentries<1000)markError(line[1],&r);
  sprintf(line[2],"%.2f +- %.2f",mean,emean);
  if(fabs(mean-m_mean)>1)markError(line[2],&r);
  sprintf(line[3],"%.2f +- %.2f",sigma,esigma);
  if(fabs(sigma-m_sigma)>2)markError(line[3],&r);
  sprintf(line[4],"%d",delay);
  r.addtableline(line,5);
  c1.SaveAs("singletower.gif");
  m_singletower->Write();

  fit(m_singletowertile,mean,emean,sigma,esigma,nentries);
  bestdelay=(int)(mean+.5);
  delay=m_acddelay-m_tkrdelay-bestdelay;
  r.linktext(line[0],"Single Tower and Tile Events", "singletowertile.gif");
  sprintf(line[1],"%d",nentries);
  if(nentries<1000)markError(line[1],&r);
  sprintf(line[2],"%.2f +- %.2f",mean,emean);
  if(fabs(mean-m_mean)>1)markError(line[2],&r);
  sprintf(line[3],"%.2f +- %.2f",sigma,esigma);
  if(fabs(sigma-m_sigma)>2)markError(line[3],&r);
  sprintf(line[4],"%d",delay);
  r.addtableline(line,5);
  c1.SaveAs("singletowertile.gif");
  m_singletowertile->Write();

  fit(m_cno,mean,emean,sigma,esigma,nentries);
  bestdelay=(int)(mean+.5);
  delay=m_acddelay-m_tkrdelay-bestdelay;
  r.linktext(line[0],"CNO", "cno.gif");
  sprintf(line[1],"%d",nentries);
  if(nentries<100)markError(line[1],&r);
  sprintf(line[2],"%.2f +- %.2f",mean,emean);
  if(fabs(mean-m_mean)>3)markError(line[2],&r);
  sprintf(line[3],"%.2f +- %.2f",sigma,esigma);
  if(fabs(sigma-m_sigma)>3)markError(line[3],&r);
  sprintf(line[4],"%d",delay);
  r.addtableline(line,5);
  c1.SaveAs("cno.gif");
  m_cno->Write();

  r.endtable();

  r.newheadline("Results by Tower");
  char* towertable[]={"Tower","Number of events","Mean +- error","Sigma +- error","Optimal delay difference"};
  r.starttable(towertable,5);
  for (int i=0;i<16;i++){
    fit(m_bytower[i],mean,emean,sigma,esigma,nentries);
    bestdelay=(int)(mean+.5);
    delay=m_acddelay-m_tkrdelay-bestdelay;
    sprintf(name,"Tower %d",i);
    sprintf(title,"tower_%d.gif",i);
    r.linktext(line[0],name,title);
    sprintf(line[1],"%d",nentries);
    if(nentries<100)markError(line[1],&r);
    sprintf(line[2],"%.2f +- %.2f",mean,emean);
    if(fabs(mean-m_mean)>1)markError(line[2],&r);
    sprintf(line[3],"%.2f +- %.2f",sigma,esigma);
    if(fabs(sigma-m_sigma)>2)markError(line[3],&r);
    sprintf(line[4],"%d",delay);
    r.addtableline(line,5);
    sprintf(name,"tower_%d.gif",i);
    c1.SaveAs(name);
    m_bytower[i]->Write();
  }
  r.endtable();
    
  r.newheadline("Results by Tile");
  char* tiletable[]={"Tile","Number of events","Mean +- error","Sigma +- error","Optimal delay difference"};
  r.starttable(tiletable,5);
  for (int i=0;i<5;i++){
    for (int j=0;j<5;j++){
      fit(m_bytile[AcdId::indexFromTile(i*10+j)],mean,emean,sigma,esigma,nentries);
      bestdelay=(int)(mean+.5);
      delay=m_acddelay-m_tkrdelay-bestdelay;
      sprintf(name,"Tile %03d",i*10+j);
      sprintf(title,"tile_%03d.gif",i*10+j);
      r.linktext(line[0],name,title);
      sprintf(line[1],"%d",nentries);
      if(nentries<50)markError(line[1],&r);
      sprintf(line[2],"%.2f +- %.2f",mean,emean);
      if(fabs(mean-m_mean)>1)markError(line[2],&r);
      sprintf(line[3],"%.2f +- %.2f",sigma,esigma);
      if(fabs(sigma-m_sigma)>2)markError(line[3],&r);
      sprintf(line[4],"%d",delay);
      r.addtableline(line,5);
      sprintf(name,"tile_%03d.gif",i*10+j);
      c1.SaveAs(name);
      m_bytile[AcdId::indexFromTile(i*10+j)]->Write();
    }
  }
  for (int l=1;l<5;l++){
    for (int i=0;i<2;i++){
      for (int j=0;j<5;j++){
	fit(m_bytile[AcdId::indexFromTile(l*100+i*10+j)],mean,emean,sigma,esigma,nentries);
	bestdelay=(int)(mean+.5);
	delay=m_acddelay-m_tkrdelay-bestdelay;
	sprintf(name,"Tile %03d",l*100+i*10+j);
	sprintf(title,"tile_%03d.gif",l*100+i*10+j);
	r.linktext(line[0],name,title);
	sprintf(line[1],"%d",nentries);
        if(nentries<50)markError(line[1],&r);
	sprintf(line[2],"%.2f +- %.2f",mean,emean);
	if(fabs(mean-m_mean)>1)markError(line[2],&r);
	sprintf(line[3],"%.2f +- %.2f",sigma,esigma);
	if(fabs(sigma-m_sigma)>2)markError(line[3],&r);
	sprintf(line[4],"%d",delay);
	r.addtableline(line,5);
	sprintf(name,"tile_%03d.gif",l*100+i*10+j);
	c1.SaveAs(name);
	m_bytile[AcdId::indexFromTile(l*100+i*10+j)]->Write();
      }
    }
  }
  r.endtable();
  r.newheadline("Status");
  if (m_status==0)r.addstatus("Passed");
  else r.addstatus("Warning");
  r.writereport();
  f.Close();
}
//
// Main event loop
//
void treqACD::Go(Long64_t numEvents)
{        
    // For this dump example, activate read in of all tree branches

    ActivateAllBranches(); 

    // To speed up read in, can choose to disable some branches 
    digiTree->SetBranchStatus("*",0);  // enable all branches
    digiTree->SetBranchStatus("m_gem",1); 
    digiTree->SetBranchStatus("m_levelOneTrigger",1); 
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
        
        GetEvent(ievent);

        if(ievent%1000==0) 
            std::cout << "** Processing Event " << ievent << std::endl;  
	const Gem gem=evt->getGem();
	UShort_t condarracd=gem.getCondArrTime().roi();
	UShort_t condarrcno=gem.getCondArrTime().cno();
	UShort_t condarrtkr=gem.getCondArrTime().tkr();
	//if (evt->getL1T().getGemEngine()==5)continue; // engine 5 is heavily prescaled
	if(gem.getConditionSummary()&0x10){ //CNO
	  m_cno->Fill((double)condarrcno-condarrtkr);
	}
	if(gem.getConditionSummary()!=3)continue;  // TKR && ROI
/*
	int miplayers[8];
	for (int i=0;i<8;i++)miplayers[i]=0;
	if (rec){
	  TObjArray *xtalRecCol = rec->getCalRecon()->getCalXtalRecCol();
	  TIter xtalIter(xtalRecCol);
	  CalXtalRecData *xtal = 0;
	  if (xtalRecCol){
	    if(xtalRecCol->GetEntries()>0 ){
	      while ((xtal = (CalXtalRecData*)xtalIter.Next())) {
		const idents::CalXtalId cId(xtal->getPackedId());
	      	Int_t tower  = cId.getTower();
	      	Int_t layer  = cId.getLayer();
		double xenergy=xtal->getEnergy();
		if (xenergy>5&&xenergy<25)miplayers[layer]=1;
		const XtalIdx xtalIdx(cId);
	      }
	    }
	  }
	}
	int nlayers=0;
	for (int i=0;i<8;i++)nlayers+=miplayers[i];
	if(nlayers<3)continue; // MIP filter
*/
	m_allevents->Fill((double)condarracd);
	std::vector<AcdId>tilelist;
	const GemTileList acdTiles = gem.getTileList();
	UInt_t acdwordxy=acdTiles.getXy();
	UInt_t acdwordxzm=acdTiles.getXzm();
	UInt_t acdwordxzp=acdTiles.getXzp();
	UInt_t acdwordyzp=acdTiles.getYzp();
	UInt_t acdwordyzm=acdTiles.getYzm();
	UInt_t acdwordrbn=acdTiles.getRbn();
	// fill the veto tile list
	for (int i=0;i<25;i++){
	  if (acdwordxy&(1<<i))tilelist.push_back(AcdId(i/5*10+i%5));
	  if (i<16){
	    if (acdwordyzm&(1<<i))tilelist.push_back(AcdId(100+i/5*10+i%5));
	    if (acdwordxzm&(1<<i))tilelist.push_back(AcdId(200+i/5*10+i%5));
	    if (acdwordyzp&(1<<i))tilelist.push_back(AcdId(300+i/5*10+i%5));
	    if (acdwordxzp&(1<<i))tilelist.push_back(AcdId(400+i/5*10+i%5));
	  }
	  if (i<8){
	    if (acdwordrbn&(1<<i))tilelist.push_back(AcdId((5+i/4)*100+i%4));
	  }
	}
	UShort_t tkrbits=gem.getTkrVector();
	std::vector<int> tkrvector;
	for (int i=0;i<16;i++){
	  if((1<<i)&tkrbits)tkrvector.push_back(i);
	}
	if (rec){
	  AcdRecon *acdRec = rec->getAcdRecon();
	  if (acdRec){
	    int nhit=tilelist.size();
	    bool intersect=false;
	    for (int i=0;i<nhit;i++){
	      for (unsigned j=0;j<acdRec->nAcdIntersections();j++){
		if (tilelist[i]==acdRec->getAcdTkrIntersection(j)->getTileId()){
		  //		  pathlength=acdRec->getAcdTkrIntersection(j)->getPathLengthInTile()/10.;
		  intersect=true;
		  break;
		}
	      }
	    }	
	    if (intersect==true){
	      m_intersect->Fill(condarracd);
	      if(tkrvector.size()==1){
		m_singletower->Fill(condarracd);
		m_bytower[tkrvector[0]]->Fill(condarracd);
	      }
	      if (tilelist.size()==1){
		m_singletile->Fill(condarracd);
		if(tkrvector.size()==1)m_singletowertile->Fill(condarracd);
		m_bytile[AcdId::indexFromTile(tilelist[0].getId())]->Fill(condarracd);
	      }
	    }
	  }
	}
    }  // end analysis code in event loop
    

std::cout << "**** End of event loop **** " << std::endl;


}

void treqACD::markError(char* line, TestReport *r){
  char text[128];
  m_status=1;
  strcpy(text,line);
  r->redtext(line,text); 
}
   
#endif
