/** 
*@class vetoeff
*
*@brief Root Tree Ananlysis Trigger time structure of a run
*
* This class is for the analysis of veto efficiency runs
*
* Version 0.1 18-Oct-2006 Martin Kocian Initial version
*/

#ifndef vetoeff_cxx
#define vetoeff_cxx 1

#include "calibGenTRG/vetoeff.h"
#include "TestReport.h"
#include "eval.h"

#if !defined(__CINT__)
// Need these includes if we wish to compile this code
#include "calibGenTRG/RootTreeAnalysis.h"
#include "TProfile.h"
#include "TCanvas.h"
#include <time.h>
#else
class RTAgeom;
class RootTreeAnalysis;
#endif

ClassImp(vetoeff)
  
// 
// Initialization
//
vetoeff::vetoeff() 
{
}

void vetoeff::inithistos(){

  char histname[128];
  char histtitle[128];
  m_baddigievents=0;
  m_tilet=new TH1D("tilet","Tiles in GEM",600,-.5,599.5);
  m_tilet->GetXaxis()->SetTitle("side*100+row*10+column");
  m_tiled=new TH1D("tiled","Tiles in Digis",600,-.5,599.5);
  m_tiled->GetXaxis()->SetTitle("side*100+row*10+column");
  for (int i=0;i<16;i++){
    m_goodcount[i]=0;
    m_badcount[i]=0;
    m_goodcountev[i]=0;
    m_baddigicount[i]=0;
    sprintf(histname,"good_ACD_tiles_tower_%d" ,i);
    sprintf(histtitle,"Good ACD tiles tower_%d", i);
    m_goodtiles[i]=new TH1D(histname,histtitle,600,-.5,599.5);
    m_goodtiles[i]->GetXaxis()->SetTitle("side*100+row*10+column");
    sprintf(histname,"bad_ACD_tiles_tower_%d" ,i);
    sprintf(histtitle,"Bad ACD tiles tower_%d", i);
    m_badtiles[i]=new TH1D(histname,histtitle,600,-.5,599.5);
    m_badtiles[i]->GetXaxis()->SetTitle("side*100+row*10+column");
    sprintf(histname,"bad_ACD_tiles_digi_tower_%d" ,i);
    sprintf(histtitle,"Bad ACD tiles in digis tower_%d", i);
    m_baddigis[i]=new TH1D(histname,histtitle,600,-.5,599.5);
    m_baddigis[i]->GetXaxis()->SetTitle("side*100+row*10+column");
  }
  roicond=new TH1D("roicond","Conditions summary with ROI set",256,-.5,255.5);
  digicond=new TH1D("digicond","Conditions summary with digi shadow",256,-.5,255.5);
  digiglt=new TH1D("digiglt","Glt word with digi shadow",256,-.5,255.5);
  glt3=new TH1D("glt3","Glt==3, cond !=3, shadow event?",2,-.5,1.5);

}  

void vetoeff::savehistos(char* filename){
  strcpy(m_filename,filename);
  TFile f(filename,"recreate");
  for (int i=0;i<16;i++){
    m_goodtiles[i]->Write();
    m_badtiles[i]->Write();
    m_baddigis[i]->Write();
  }
  roicond->Write();
  digicond->Write();
  digiglt->Write();
  glt3->Write();
  m_tilet->Write();
  m_tiled->Write();
  f.Close();
}
  

void vetoeff::readParameterFile(char* pf){
  ifstream ifs(pf);
  int num;
  char c;
  int i=0;
  while (!ifs.eof()){
    do{
      c=ifs.get();
      if (c=='\n')i++;
    }while(c==' ');
    ifs.unget();
    ifs>>num;
    m_roi[i].push_back(num);
  }
}


void vetoeff::writeoutresults(char* report){
  std::cout<<"Gem shadow and also digi shadow: "<<m_gemdigi<<std::endl;
  std::cout<<"Gem shadow but no digi shadow: "<<m_gemnodigi<<std::endl;
  std::cout<<"Digi shadow but no gem shadow: "<<m_diginogem<<std::endl;
  float digieff=efficiency(m_gemdigi,m_gemdigi+m_gemnodigi);
  float edigieff=efferror(m_gemdigi,m_gemdigi+m_gemnodigi);
  float gemeff=efficiency(m_gemdigi,m_gemdigi+m_diginogem);
  float egemeff=efferror(m_gemdigi,m_gemdigi+m_diginogem);
  std::cout<<"Digi shadow efficiency "<<digieff<<" +- "<<edigieff<<std::endl;
  std::cout<<"Gem shadow efficiency "<<gemeff<<" +- "<<egemeff<<std::endl;
  char name[128];
  char title[128];
  char* line[7];
  for (int j=0;j<7;j++)line[j]=new char[64];
  TestReport r(report);
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
  r.additem("Note","Click on links to see plots");

  gROOT->SetStyle("Plain");
/*
  TCanvas c1;
  roicond->Draw();
  c1.SaveAs("condsummary.gif");
*/
  sprintf(name,"%d",m_roicount);
  r.linktext(title,name,"condsummary.gif");
  r.additem("Triggers with ROI set",title);
  sprintf(name,"%d",m_roicount1);
  if (m_roicount1>0)r.redtext(title,name);
  else r.greentext(title,name);
  r.additem("Triggers with ROI only",title);
  sprintf(name,"%d",m_roicount3);
  if (m_roicount3>0)r.redtext(title,name);
  else r.greentext(title,name);
  r.additem("Triggers with ROI and TKR only",title);
  sprintf(name,"%d",m_baddigievents);
  if (m_baddigievents>0)r.redtext(title,name);
  else r.greentext(title,name);
  r.additem("Events with shadow in digis",title);
  sprintf(name,"%d",m_roicountveto);
  if (m_roicountveto>0)r.redtext(title,name);
  else r.greentext(title,name);
  r.additem("Events with throttle in digis",title);
  sprintf(name,"%d",m_roicount32);
  if (m_roicount32>0)r.redtext(title,name);
  else r.greentext(title,name);
  r.additem("Events with ROI only in digis",title);
  sprintf(name,"%d",m_roicount35);
  if (m_roicount35>0)r.redtext(title,name);
  else r.greentext(title,name);
  r.additem("Events with ROI and TKR in digis",title);
  sprintf(name,"%d",m_32shadow);
  r.additem("Events with GltWord==3 && GemConditionWord==2 and shadow",name);
  sprintf(name,"%d",m_32nofilter);
  r.additem("Events with GltWord==3 && GemConditionWord==2 and no filter tile",name);
  sprintf(name,"%d",m_32tkrintersec);
  r.additem("Events with GltWord==3 && GemConditionWord==2 and track intersects tile",name);
  sprintf(name,"%d",m_32noise);
  r.additem("Events with GltWord==3 && GemConditionWord==2 and random tile",name);
  char* restable[]={"Tower","Events with non-shadowing tiles","Number of non-shadowing tiles","Number of shadowing tiles"};
  r.starttable(restable,4);
  for (int i=0;i<16;i++){
    sprintf(title,"%d",m_goodcount[i]);
    if(m_goodcount[i]!=0){
      //     m_goodtiles[i]->Draw();
      sprintf(name,"goodtiles_tower_%d.gif",i);
            //c1.SaveAs(name);
      r.linktext(line[2],title,name);
    } else sprintf(line[2],title);
    sprintf(name,"%d",m_badcount[i]);
    if (m_badcount[i]>0)r.redtext(title,name);
    else r.greentext(title,name);
    if(m_badcount[i]!=0){
            //m_badtiles[i]->Draw();
      sprintf(name,"badtiles_tower_%d.gif",i);
            //c1.SaveAs(name);
      r.linktext(line[3],title,name);
    }else sprintf(line[3],title);  
    sprintf(line[0],"%d",i);
    sprintf(line[1],"%d",m_goodcountev[i]);
    r.addtableline(line,4);
  }    
  r.endtable();
//  r.newheadline("Plots for shadowing tiles");
//  for (int i=0;i<16;i++){
//    if(m_badcount[i]!=0){
//      m_badtiles[i]->Draw();
}
// Main event loop
//
void vetoeff::Go(Long64_t numEvents)
{        
    // For this dump example, activate read in of all tree branches

    ActivateAllBranches(); 

    // To speed up read in, can choose to disable some branches 
    digiTree->SetBranchStatus("*",1);  // enable all branches
    reconTree->SetBranchStatus("*",1);  // disable recon branches
    
       
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
    m_roicount=m_roicount1=m_roicount3=m_roicount32=m_roicount35=m_roicountveto=0;
    m_32shadow=m_32nofilter=m_32tkrintersec=m_32noise=0;
    m_diginogem=m_gemdigi=m_gemnodigi=0;
    m_status=0;
    
    GetEvent(m_StartEvent);
    m_gid=evt->getMetaEvent().run().id();
    std::cout<<"gid "<<m_gid<<std::endl;
    m_nev=nEvtMax-m_StartEvent;

    //}
    // BEGINNING OF EVENT LOOP
    std::vector<int>tilelist;
    std::vector<int>acddigilist;
    for (Long64_t ievent=m_StartEvent; ievent<nEvtMax; ievent++) {
        
        if (mc) mc->Clear();
        if (evt) evt->Clear();
        if (rec) rec->Clear();
	tilelist.clear();
	acddigilist.clear();
        
        GetEvent(ievent);

        if(ievent%1000==0) 
            std::cout << "** Processing Event " << ievent << std::endl;  
	const Gem GemData  = evt->getGem();
	if (GemData.getDeltaEventTime()<2000)continue;
	const GemTileList acdTiles = GemData.getTileList();
	UInt_t acdwordxy=acdTiles.getXy();
	UInt_t acdwordxzm=acdTiles.getXzm();
	UInt_t acdwordxzp=acdTiles.getXzp();
	UInt_t acdwordyzp=acdTiles.getYzp();
	UInt_t acdwordyzm=acdTiles.getYzm();
	UInt_t acdwordrbn=acdTiles.getRbn();
	 acdwordrbn=0;
	// fill the veto tile list
	for (int i=0;i<25;i++){
	  if (acdwordxy&(1<<i))tilelist.push_back(i/5*10+i%5);
	  if (i<16){
	    if (acdwordyzm&(1<<i))tilelist.push_back(100+i/5*10+i%5);
	    if (acdwordxzm&(1<<i))tilelist.push_back(200+i/5*10+i%5);
	    if (acdwordyzp&(1<<i))tilelist.push_back(300+i/5*10+i%5);
	    if (acdwordxzp&(1<<i))tilelist.push_back(400+i/5*10+i%5);
	  }
	  if (i<8){
	    if (acdwordrbn&(1<<i))tilelist.push_back((5+i/4)*100+i%4);
	  }
	}
	AcdRecon *acdRec = rec->getAcdRecon();
	if (acdRec){
	  int nhit=acdRec->nAcdHit();
	  for (int i=0;i<nhit;i++){
	    const AcdHit* acdHit=acdRec->getAcdHit(i);
	    double acdenergy=(double)acdHit->getMips();
            if (acdenergy<m_thresh)continue;
	    AcdId id=acdHit->getId();
	    if (id.isTile() && !id.getNa())acddigilist.push_back(id.getFace()*100+id.getRow()*10+id.getColumn());
	  }	
	}
	//	const TObjArray* acdDigiCol = evt->getAcdDigiCol();
		//if (acdDigiCol){
		//  TIter acdIter(acdDigiCol);
		//  AcdDigi* acddigi=0;
		//  char tilename[10];
		//  while ((acddigi = (AcdDigi*)acdIter.Next())) {
		//    strcpy(tilename,acddigi->getTileName());
		//    if (strstr(tilename,"NA")==NULL)acddigilist.push_back(atoi(tilename));
		//  }
		//}
	UShort_t cond=GemData.getConditionSummary();
	if (cond&0x20)continue;
	UInt_t gltword=evt->getL1T( ).getTriggerWord( );
	gltword &=0xff;
	if ((gltword==3 ||gltword==2)&& cond==2){
	  for (int i=0;i<acddigilist.size();i++){
	    m_tiled->Fill(acddigilist[i]);
	  }
	  for (int i=0;i<tilelist.size();i++){
	    m_tilet->Fill(tilelist[i]);
	  }
	}
	    
	if (cond&1){
	  roicond->Fill(cond);
	  m_roicount++;
	  if (cond==1)m_roicount1++;
	  if (cond==3)m_roicount3++;
	}
	if (evt->getL1T().getTriggerWord()&32)m_roicountveto++;
	if (evt->getL1T().getTriggerWord()==32)m_roicount32++;
	if (evt->getL1T().getTriggerWord()==35)m_roicount35++;

	UShort_t tkrvector=GemData.getTkrVector();
	if (!(cond&1)){
	  for (int i=0;i<16;i++){
	    int goodev=0;
	    if (tkrvector&(1<<i)){
	      for (int j=0;j<tilelist.size();j++){
		bool foundtile=false;
		for (int k=0;k<m_roi[i].size();k++){
		  if(tilelist[j]==m_roi[i][k]){
		    foundtile=true;
		    break;
		  }
		}
		if (foundtile==true){
		  m_badtiles[i]->Fill(tilelist[j]);
		  m_badcount[i]++;
		  m_status=1;
		}else{
		  m_goodtiles[i]->Fill(tilelist[j]);
		  m_goodcount[i]++;
		  goodev=1;
		}
	      }
	    }
	    m_goodcountev[i]+=goodev;
	  }
	}
 
	int badev=0;
	if (!(gltword&32)){
	  for (int i=0;i<16;i++){
	    UShort_t tribit=evt->getL1T( ).getDigiTriRowBits(i);
	    if (tribit){
	      //std::cout<<"trirowbits"<<std::endl;
	      for (int j=0;j<acddigilist.size();j++){
		bool foundtile=false;
		for (int k=0;k<m_roi[i].size();k++){
		  if(acddigilist[j]==m_roi[i][k]){
		    foundtile=true;
		    break;
		  }
		}
		if (foundtile==true){
		  if (cond<4)digiglt->Fill(gltword);
		  digicond->Fill(cond);
		  //std::cout<<"bad event "<<ievent<<std::endl;
		  //std::cout<<"Gem Word "<<cond<<std::endl;
		  //std::cout<<"Glt Word "<<gltword<<std::endl;
		  //std::cout<<"tower "<<i<<std::endl;
		  //std::cout<<"tile "<<acddigilist[j]<<std::endl;
		  m_baddigis[i]->Fill(acddigilist[j]);
		  m_baddigicount[i]++;
		  badev=1;
		  //m_status=1;
		}
	      }
	    }
	  }
	  m_baddigievents+=badev;
	}
	int filtertile=0;
	for (int j=0;j<acddigilist.size();j++){
	  if (acddigilist[j]<100 || acddigilist[j]%100<15)filtertile=1;
	}
	int tkrintersec=0;
	if (acdRec){
	  int nint=acdRec->nAcdIntersections();
	  for (int i=0;i<nint;i++){
	    const AcdTkrIntersection* acdint=acdRec->getAcdTkrIntersection(i);
	    AcdId id=acdint->getTileId();
	    if (id.isTile() && !id.getNa()){
	      int exptile=id.getFace()*100+id.getRow()*10+id.getColumn();
	      if(std::find(acddigilist.begin(),acddigilist.end(),exptile)!=acddigilist.end())tkrintersec=1;
	    }
	  }
	}
	    
	if(cond&1){
	  if (badev)m_gemdigi++;
	  else m_gemnodigi++;
	}
	if(badev){
	  if (!(cond&1))m_diginogem++;
	}
	if (gltword==3 && (cond&3)!=3){
	  //glt3->Fill(badev);
	  //std::cout<<"event "<<ievent<<std::endl;
	  if (badev){
	   // std::cout<<"shadow event"<<std::endl;
	    m_32shadow++;
	  }
	  else if (!filtertile){
	   // std::cout<<"no filter tile hit event"<<std::endl;
	    m_32nofilter++;
	  }
	  else if (tkrintersec){
	    //std::cout<<"tkr with tile intersection in event "<<ievent<<" "<<cond<<" "<<GemData.getRoiVector()<<std::endl;
	    m_32tkrintersec++;
	  }
	  else {
	   // std::cout<<"noise event"<<std::endl;
	    m_32noise++;
	  }
	}
		  
  
	  
	
	
	
    }  // end analysis code in event loop
    

std::cout << "**** End of event loop **** " << std::endl;

// GEM time counter summary
//igemUtil->dumpTimeCounter();

// Endjob histograms  
//AcdHend();

}

float vetoeff::efficiency(int a,int b){
  if (b>0){
    return (float)a/(float)b;
  } else {
    return 0;
  }
}
float vetoeff::efferror(int a,int b){
  float eff;
  if (b>0){
    eff= (float)a/(float)b;
    return sqrt(eff*(1.-eff)/(float)b);
  } else {
    return 0;
  }
}



#endif
