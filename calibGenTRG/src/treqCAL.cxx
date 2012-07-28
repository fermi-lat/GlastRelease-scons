/** 
*@class takedata_mt
*
*@brief Root Tree Ananlysis Trigger time structure of a run
*
* This class is for Trigger Root Tree Ananlysis of the run time structure
*
* Version 0.1 30-Jun-2006 Martin Kocian Initial version
*/


#include "calibGenTRG/treqCAL.h"
#include "TestReport.h"
#include "eval.h"
#include "tackcuts.h"

#include <TF1.h>
#include <TCanvas.h>
#include "calibGenTRG/RootTreeAnalysis.h"
#include "TProfile.h"
#include "TMath.h"

#include "CalGeom.h"

#include <list>
#include <time.h>

using namespace CLHEP;
using namespace CalUtil;

ClassImp(treqCAL)
  
// 
// Initialization
//
treqCAL::treqCAL() : m_caltuplefile(0),m_caltuple(0), m_useKalman(false), m_useToT(false), m_useMip(false), m_useGamma(false),
		     m_kalmanLower(700), m_ToTLower(1), m_ToTUpper(1.6), m_MipLower(3)
{
}

void treqCAL::initCalTuple(const char* filename){
  if (m_caltuplefile){
    delete m_caltuple;
    delete m_caltuplefile;
  }
    m_caltuplefile=TFile::Open(filename, "READ");
    if(m_caltuplefile->IsOpen()){
      std::cout<<"        CAL:    "<<filename<<std::endl;
      m_caltuple=(TTree*)gDirectory->Get("CalTuple");
      m_caltuple->SetBranchAddress("CalXtalFaceSignal[16][8][12][2]",m_calxtalfacesignal);
    }else{
      std::cout<<"CAL file "<<filename<<" not found"<<std::endl;
    }
}

 
void treqCAL::setParameters(int tkrdelay, int caldelay){
  m_tkrdelay=tkrdelay;
  m_caldelay=caldelay;
  std::cout<<"TKR TREQ delay: "<<m_tkrdelay<<std::endl;
  std::cout<<"CAL TREQ delay: "<<m_caldelay<<std::endl;
}

void treqCAL::inithistos(){

  char histname[128];
  char histtitle[128];
    //sprintf(histname,"ACD_adchist");
    //sprintf(histtitle,"ACD ADC histogram");
    //m_acdadchist=new TH1D(histname,histtitle,100,0,3);
    //m_acdadchist->GetXaxis()->SetTitle("MIPS");
  m_clallevents=new TH1D("clallevents","All Events with CAL LO",63,-31.5,31.5);
  m_clallevents->GetXaxis()->SetTitle("CAL LO conditions arrival time (ticks)");
  m_clvsenergy=new TH2D("clvsenergy","CAL LO arrival time vs. highest crystal energy",100,0,2000,63,-31.5,31.5);
  m_clvsenergy->GetYaxis()->SetTitle("CAL LO conditions arrival time (ticks)");
  m_clvsenergy->GetXaxis()->SetTitle("Energy (MeV)");
  m_chvsenergy=new TH2D("chvsenergy","CAL HI arrival time vs. highest crystal energy",100,0,10000,63,-31.5,31.5);
  m_chvsenergy->GetYaxis()->SetTitle("CAL HI conditions arrival time (ticks)");
  m_chvsenergy->GetXaxis()->SetTitle("Energy (MeV)");
  m_challevents=new TH1D("challevents","All Events with CAL HI",63,-31.5,31.5);
  m_challevents->GetXaxis()->SetTitle("CAL HI conditions arrival time (ticks)");
  m_clintersect=new TH1D("clintersect","Events where Track intersects highest energy crystal" ,63,-31.5,31.5);
  m_clintersect->GetXaxis()->SetTitle("CAL LO conditions arrival time (ticks)");
  m_chintersect=new TH1D("chintersect","Events where Track intersects highest energy crystal" ,63,-31.5,31.5);
  m_chintersect->GetXaxis()->SetTitle("CAL HI conditions arrival time (ticks)");
  m_clsingletower=new TH1D("clsingletower","Single Tower Events",63,-31.5,31.5);
  m_clsingletower->GetXaxis()->SetTitle("CAL LO conditions arrival time (ticks)");
  m_chsingletower=new TH1D("chsingletower","Single Tower Events",63,-31.5,31.5);
  m_chsingletower->GetXaxis()->SetTitle("CAL HI conditions arrival time (ticks)");
  for (int i=0;i<16;i++){
    sprintf(histname,"cltower_%d",i);
    sprintf(histtitle,"Tower %d CAL LO",i);
    m_clbytower[i]=new TH1D(histname,histtitle,63,-31.5,31.5);
    m_clbytower[i]->GetXaxis()->SetTitle("CAL LO conditions arrival time (ticks)");
  }
  for (int i=0;i<16;i++){
    sprintf(histname,"chtower_%d",i);
    sprintf(histtitle,"Tower %d CAL HI",i);
    m_chbytower[i]=new TH1D(histname,histtitle,63,-31.5,31.5);
    m_chbytower[i]->GetXaxis()->SetTitle("CAL HI conditions arrival time (ticks)");
  }
}  


void treqCAL::fit(TH1D* histo,double& mean, double &emean, double& sigma, double& esigma, int& nevents){
  histo->Fit("gaus");
  mean=histo->GetFunction("gaus")->GetParameter(1);
  emean=histo->GetFunction("gaus")->GetParError(1);
  sigma=histo->GetFunction("gaus")->GetParameter(2);
  esigma=histo->GetFunction("gaus")->GetParError(2);
  nevents=(int)histo->GetEntries();
}

void treqCAL::getMean(TH1D* histo,double& mean, double &emean, double& sigma, double& esigma, int& nevents){
  nevents=(int)histo->GetEntries();
  mean=histo->GetMean();
  emean=histo->GetMeanError();
  sigma=histo->GetRMS();
  esigma=histo->GetRMSError();
  histo->Draw();
}
  

void treqCAL::writeoutresults(const char* reportname, const char* filename){
  m_status=0;
  char name[128];
  char title[128];
  char* line[7];
  for (int j=0;j<7;j++)line[j]=new char[512];
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(0);
  gStyle->SetPalette(1);
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
  r.additem("Application", "treqCAL");
  r.additem("Description", "Analysis of a Treq_TkrCal run");
  r.additem("Type of analysis","Offline");
  sprintf (name,"%d",m_gid);
  r.additem("Run ground id",name);
  sprintf (name,"%d",m_nev);
  r.additem("Number of events",name);
  sprintf (name,"%d",m_tkrdelay);
  r.additem("TKR delay",name);
  sprintf (name,"%d",m_caldelay);
  r.additem("ACD delay",name);

  r.newheadline("CALLO General Results");
  r.addimage("clallevents.gif");
  m_clvsenergy->Draw("colz");
  c1.SaveAs("clvsenergy.gif");
  r.addimage("clvsenergy.gif");
  m_clvsenergy->Write();
  m_clvsenergy->FitSlicesY();
  TH1D* clslices=(TH1D*)gDirectory->Get("clvsenergy_1");
  clslices->SetAxisRange(-15,15,"y");
  clslices->Draw();
  c1.SaveAs("clslices.gif");
  r.addimage("clslices.gif");
  clslices->Write();

  char* restable[]={"Cut","Number of events","Mean +- error","Sigma +- error","Optimal delay difference"};
  double mean,emean,sigma,esigma;
  int nentries, bestdelay, delay;
  r.starttable(restable,5);
  getMean(m_clallevents,mean,emean,sigma,esigma,nentries);
  m_clmean=mean;
  m_clsigma=sigma;
  bestdelay=TMath::Nint(mean);
  delay=m_tkrdelay-m_caldelay+bestdelay;
  r.linktext(line[0],"All Events (CAL LO)", "clallevents.gif");
  sprintf(line[1],"%d",nentries);
  if(nentries<1000)markError(line[1],&r);
  sprintf(line[2],"%.2f +- %.2f",mean,emean);
  sprintf(line[3],"%.2f +- %.2f",sigma,esigma);
  sprintf(line[4],"%d",delay);
  r.addtableline(line,5);
  c1.SaveAs("clallevents.gif");
  m_clallevents->Write();

  getMean(m_clintersect,mean,emean,sigma,esigma,nentries);
  bestdelay=TMath::Nint(mean);
  delay=m_tkrdelay-m_caldelay+bestdelay;
  r.linktext(line[0],"Events with TKR Max E Crystal Intersection (CAL LO)", "clintersect.gif");
  sprintf(line[1],"%d",nentries);
  if(nentries<1000)markError(line[1],&r);
  sprintf(line[2],"%.2f +- %.2f",mean,emean);
  if(fabs(mean-m_clmean)>1)markError(line[2],&r);
  sprintf(line[3],"%.2f +- %.2f",sigma,esigma);
  if(fabs(sigma-m_clsigma)>2)markError(line[3],&r);
  sprintf(line[4],"%d",delay);
  r.addtableline(line,5);
  c1.SaveAs("clintersect.gif");
  m_clintersect->Write();

  getMean(m_clsingletower,mean,emean,sigma,esigma,nentries);
  bestdelay=TMath::Nint(mean);
  delay=m_tkrdelay-m_caldelay+bestdelay;
  r.linktext(line[0],"Single Tower Events (CAL LO)", "clsingletower.gif");
  sprintf(line[1],"%d",nentries);
  if(nentries<1000)markError(line[1],&r);
  sprintf(line[2],"%.2f +- %.2f",mean,emean);
  if(fabs(mean-m_clmean)>1)markError(line[2],&r);
  sprintf(line[3],"%.2f +- %.2f",sigma,esigma);
  if(fabs(sigma-m_clsigma)>2)markError(line[3],&r);
  sprintf(line[4],"%d",delay);
  r.addtableline(line,5);
  c1.SaveAs("clsingletower.gif");
  m_clsingletower->Write();

  r.endtable();

  r.newheadline("CALLO Results by Tower");
  char* towertable[]={"Tower","Number of events","Mean +- error","Sigma +- error","Optimal delay difference"};
  r.starttable(towertable,5);
  for (int i=0;i<16;i++){
    getMean(m_clbytower[i],mean,emean,sigma,esigma,nentries);
    bestdelay=TMath::Nint(mean);
    delay=m_tkrdelay-m_caldelay+bestdelay;
    sprintf(name,"Tower %d (CAL LO)",i);
    sprintf(title,"cltower_%d.gif",i);
    r.linktext(line[0],name,title);
    sprintf(line[1],"%d",nentries);
    if(nentries<100)markError(line[1],&r);
    sprintf(line[2],"%.2f +- %.2f",mean,emean);
    if(fabs(mean-m_clmean)>1)markError(line[2],&r);
    sprintf(line[3],"%.2f +- %.2f",sigma,esigma);
    if(fabs(sigma-m_clsigma)>2)markError(line[3],&r);
    sprintf(line[4],"%d",delay);
    r.addtableline(line,5);
    sprintf(name,"cltower_%d.gif",i);
    c1.SaveAs(name);
    m_clbytower[i]->Write();
  }
  r.endtable();

 
  r.newheadline("CALHI General Results");
  r.addimage("challevents.gif");
  m_chvsenergy->Draw("colz");
  c1.SaveAs("chvsenergy.gif");
  r.addimage("chvsenergy.gif");
  m_chvsenergy->Write();
  m_chvsenergy->FitSlicesY();
  TH1D* chslices=(TH1D*)gDirectory->Get("chvsenergy_1");
  chslices->SetAxisRange(-15,15,"y"); 
  chslices->Draw();
  c1.SaveAs("chslices.gif");
  r.addimage("chslices.gif");
  chslices->Write();

  r.starttable(restable,5);
  
  getMean(m_challevents,mean,emean,sigma,esigma,nentries);
  m_chmean=mean;
  m_chsigma=sigma;
  bestdelay=TMath::Nint(mean);
  delay=m_tkrdelay-m_caldelay+bestdelay;
  r.linktext(line[0],"All Events (CAL HI)", "challevents.gif");
  sprintf(line[1],"%d",nentries);
  if(nentries<1000)markError(line[1],&r);
  sprintf(line[2],"%.2f +- %.2f",mean,emean);
  sprintf(line[3],"%.2f +- %.2f",sigma,esigma);
  sprintf(line[4],"%d",delay);
  r.addtableline(line,5);
  c1.SaveAs("challevents.gif");
  m_challevents->Write();

  getMean(m_chintersect,mean,emean,sigma,esigma,nentries);
  bestdelay=TMath::Nint(mean);
  delay=m_tkrdelay-m_caldelay+bestdelay;
  r.linktext(line[0],"Events with TKR Max E Crystal Intersection (CAL HI)", "chintersect.gif");
  sprintf(line[1],"%d",nentries);
  if(nentries<1000)markError(line[1],&r);
  sprintf(line[2],"%.2f +- %.2f",mean,emean);
  if(fabs(mean-m_chmean)>1)markError(line[2],&r);
  sprintf(line[3],"%.2f +- %.2f",sigma,esigma);
  if(fabs(sigma-m_chsigma)>2)markError(line[3],&r);
  sprintf(line[4],"%d",delay);
  r.addtableline(line,5);
  c1.SaveAs("chintersect.gif");
  m_chintersect->Write();

  getMean(m_chsingletower,mean,emean,sigma,esigma,nentries);
  bestdelay=TMath::Nint(mean);
  delay=m_tkrdelay-m_caldelay+bestdelay;
  r.linktext(line[0],"Single Tower Events (CAL HI)", "chsingletower.gif");
  sprintf(line[1],"%d",nentries);
  if(nentries<1000)markError(line[1],&r);
  sprintf(line[2],"%.2f +- %.2f",mean,emean);
  if(fabs(mean-m_chmean)>1)markError(line[2],&r);
  sprintf(line[3],"%.2f +- %.2f",sigma,esigma);
  if(fabs(sigma-m_chsigma)>2)markError(line[3],&r);
  sprintf(line[4],"%d",delay);
  r.addtableline(line,5);
  c1.SaveAs("chsingletower.gif");
  m_chsingletower->Write();

  r.endtable();

  r.newheadline("CALHI Results by Tower");
  r.starttable(towertable,5);
  for (int i=0;i<16;i++){
    getMean(m_chbytower[i],mean,emean,sigma,esigma,nentries);
    bestdelay=TMath::Nint(mean);
    delay=m_tkrdelay-m_caldelay+bestdelay;
    sprintf(name,"Tower %d (CAL HI)",i);
    sprintf(title,"chtower_%d.gif",i);
    r.linktext(line[0],name,title);
    sprintf(line[1],"%d",nentries);
    if(nentries<100)markError(line[1],&r);
    sprintf(line[2],"%.2f +- %.2f",mean,emean);
    if(fabs(mean-m_chmean)>1)markError(line[2],&r);
    sprintf(line[3],"%.2f +- %.2f",sigma,esigma);
    if(fabs(sigma-m_chsigma)>2)markError(line[3],&r);
    sprintf(line[4],"%d",delay);
    r.addtableline(line,5);
    sprintf(name,"chtower_%d.gif",i);
    c1.SaveAs(name);
    m_chbytower[i]->Write();
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
void treqCAL::Go(Long64_t numEvents)
{        
    // For this dump example, activate read in of all tree branches

    ActivateAllBranches(); 

    // To speed up read in, can choose to disable some branches 
    digiTree->SetBranchStatus("*",0);  // enable all branches
    digiTree->SetBranchStatus("m_gem",1);  // enable all branches
    digiTree->SetBranchStatus("m_metaEvent",1);  // enable all branches
    digiTree->SetBranchStatus("m_levelOneTrigger",1);  // enable all branches
    digiTree->SetBranchStatus("m_obfFilterStatus",1);  // enable all branches
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
    int eng5=0;
    int gamma=0;

    //}
    // BEGINNING OF EVENT LOOP
    for (Long64_t ievent=m_StartEvent; ievent<nEvtMax; ievent++) {
        
        if (mc) mc->Clear();
        if (evt) evt->Clear();
        if (rec) rec->Clear();
        
        GetEvent(ievent);
	if (m_caltuple)m_caltuple->GetEvent(ievent);
        if(ievent%10000==0) 
            std::cout << "\r** Processing Event " << ievent << std::flush;  
	const Gem gem=evt->getGem();
	// require TKR && CAL but No CNO
	if( (!(gem.getConditionSummary()&0x2) || !(gem.getConditionSummary()&0xc)) && !(gem.getConditionSummary()&0x10) )continue;  
	double maxeface=0;
	bool maxveto=false;
	int miplayers[8];
	double Elayers[8];
	for (int i=0;i<8;i++){
	  miplayers[i]=0;
	  Elayers[i]  =0;
	  }
	XtalIdx maxIdxface;
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
	      	Int_t column  = cId.getColumn();
		
		double xenergy=xtal->getEnergy();
		if (xenergy>5&&xenergy<25) miplayers[layer]=1;  // Require a mip in a single xtal
		Elayers[layer]+=xenergy; // and calculate Elayer for the whole LAT
		
		const XtalIdx xtalIdx(cId);
		bool veto=false;
		if(m_calxtalfacesignal[tower][layer][column][0]>0 && m_calxtalfacesignal[tower][layer][column][1]>0){
		  float faceratio=m_calxtalfacesignal[tower][layer][column][0]/m_calxtalfacesignal[tower][layer][column][1];
		  if(faceratio<1)faceratio=1/faceratio;
		  if(faceratio>2){
		    veto=true;
		  }
		}else{
		  veto=true;
		}
		if(m_calxtalfacesignal[tower][layer][column][0]>maxeface){
		  maxeface=m_calxtalfacesignal[tower][layer][column][0];
		  maxIdxface=xtalIdx;
		  maxveto=veto;
		}
		if(m_calxtalfacesignal[tower][layer][column][1]>maxeface){
		  maxeface=m_calxtalfacesignal[tower][layer][column][1];
		  maxIdxface=xtalIdx;
		  maxveto=veto;
		}
	      } //end of while	      
	    }
	  }
	} // end of if rec ...but we're using recon for the tkr afterwards.
	if (maxeface<90 || maxveto==true)continue;
	
	// Check Layers Energy to see whether there is a mip in there or not
	// layer counts for all CAL towers, this in someway cleans more events than checking only one tower.
	for (int i=0;i<8;i++)
	   miplayers[i] = (Elayers[i]>25)?0:1;

	int nlayers=0;
	for (int i=0;i<8;i++)nlayers+=miplayers[i];
	if(m_useMip && nlayers<m_MipLower)continue; // MIP filter
	// extract track information from recon file
	TkrRecon* tkrRecon = rec->getTkrRecon();
	TObjArray* tracks = tkrRecon->getTrackCol();
	if(tracks->GetEntries()>0){
	  // Select here event that have more than 1 track
	  //if( (tracks->GetEntries())!=1 && m_useOneTrack)continue;
	  //first track
	  TkrTrack* tkrTrack = dynamic_cast<TkrTrack*>(tracks->At(0));
	  if(tkrTrack) {	    	
	    // event cut on Kalman energy 
	    if(m_useKalman && tkrTrack->getKalEnergy( )<m_kalmanLower)continue;
	    // event cut track chi square 
	    //if( (tkrTrack->getChiSquareSmooth()>20) && m_useOneTrack)continue; 
	    // and now check the ToT cut
	    if (m_useToT){
	      // calculate average TOT to cut on
	      double Tkr_1_ToTAve=0;	      
	      int clustersOnTrack=0;
	      TIterator* hitIter = tkrTrack->Iterator();
	      TObject*   hitObj  = 0;
	      while ((hitObj = hitIter->Next())!=0){
		TkrTrackHit*        hit = dynamic_cast<TkrTrackHit*>(hitObj);
		unsigned int bits = hit->getStatusBits();
		if ((bits & TkrTrackHit::HITISSSD)==0) continue;
		const TkrCluster* cluster = hit->getClusterPtr();	       
		int size =  (int) (const_cast<TkrCluster*>(cluster))->getSize();
		// get the local slopes
		double slope  = fabs(hit->getMeasuredSlope(TkrTrackHit::SMOOTHED));
		double slope1 = fabs(hit->getNonMeasuredSlope(TkrTrackHit::SMOOTHED));
		
		// theta1 is the projected angle across the strip
		double theta1       = atan(slope1);
		
		double aspectRatio = 0.228/0.400;
		double totMax      =  250.;   // counts
		double threshold   =  0.25;   // Mips
		double countThreshold = 15; // counts
		double normFactor  =  1./53.;
		
		double mips = cluster->getMips();
		if(mips<0.0 || mips>10.0) continue;
		clustersOnTrack++;
		double tot = cluster->getToT();
		if(tot>=totMax) tot = totMax;
		double path1 = 1.0;
		
		// get the path length for the hit
		// tries to get the average
		// the calculation is part analytic, part approximation and part fudge.
		//   more work is definitely in order!
		
		// theta1 first
		if (tot>=totMax) { tot = normFactor*(totMax+countThreshold); }
		else {
		  double costh1 = cos(theta1);
		  if (size==1) {
		    double sinth1 = sin(theta1);
		    if (slope1< aspectRatio) {
		      path1 = (1./costh1*(aspectRatio-slope1) + 
			       (1/costh1 - 0.5*threshold)*(2*threshold*sinth1))
			/(aspectRatio - slope1 + 2*threshold*sinth1);
		    } else if (slope1<aspectRatio/(1-2.*threshold*costh1)) {
		      path1 = 1; //1/costh1 - threshold*costh1;
		    } else { 
		      path1 = 1;
		    }
		  }
		  else if (size==2) {
		    if (slope1<aspectRatio/(1.-threshold*costh1)) {
		      path1 = 0.75/costh1 -0.5*threshold;
		    } else if (slope1<2.*aspectRatio/(1.-2*threshold*costh1)) { 
		      path1 = aspectRatio/sin(theta1);
		    } else {
		      path1 = 1.0;
		    }
		  } else {
		    if(slope1>aspectRatio/(1.- 2.*threshold*costh1)) {
		      path1 = aspectRatio/sin(theta1);
		    } else {
		      path1 = 1.0;
		    }
		  }
		  double factor = path1*costh1*slope;
		  double path2 = sqrt(path1*path1 + factor*factor);
		  mips /= path2;
		} // end of else totMax
		
		Tkr_1_ToTAve += mips;
	      }// end of while on hits
	      
	      if(clustersOnTrack>3) {
		Tkr_1_ToTAve /= clustersOnTrack;
	      }
	      if (Tkr_1_ToTAve>m_ToTUpper || Tkr_1_ToTAve<m_ToTLower)continue; // only use ToT between 1 and 1.6
	      
	    } //end of useToT
	  }  // end of if tkrTrack
	} // end of if tkrTrack GetEntries
        const LpaGammaFilter *gamma = evt->getGammaFilter();
        UInt_t state=2;
        if (gamma){
        state  = gamma->getState();
        }
        if (m_useGamma && state!=0) continue; // Gamma Filter
 
	gamma++;
	
	if (evt->getL1T().getGemEngine()==5)eng5++;
	UShort_t condarrtkr=gem.getCondArrTime().tkr();
	UShort_t condarrcalhi=gem.getCondArrTime().calHE();
	UShort_t condarrcallo=gem.getCondArrTime().calLE();

	UShort_t tkrbits=gem.getTkrVector();
	std::vector<int> tkrvector;
	for (int i=0;i<16;i++){
	  if((1<<i)&tkrbits)tkrvector.push_back(i);
	}
	UShort_t callobits=gem.getCalLeVector();
	std::vector<int> callowvector;
	for (int i=0;i<16;i++){
	  if((1<<i)&callobits)callowvector.push_back(i);
	}
	UShort_t calhibits=gem.getCalHeVector();
	std::vector<int> calhivector;
	for (int i=0;i<16;i++){
	  if((1<<i)&calhibits)calhivector.push_back(i);
	}
	if(gem.getConditionSummary()&0x4 && maxeface>90){
	  m_clallevents->Fill((double)condarrcallo-condarrtkr);
	  m_clvsenergy->Fill(maxeface,(double)condarrcallo-condarrtkr);
	}
	if(gem.getConditionSummary()&0x8 && maxeface>900){
	  m_challevents->Fill((double)condarrcalhi-condarrtkr);
	  m_chvsenergy->Fill(maxeface,(double)condarrcalhi-condarrtkr);
	}
	if (tkrvector.size()==1){
	  if(callowvector.size()==1 && tkrvector[0]==callowvector[0] && maxIdxface.getTwr()==callowvector[0] && maxeface>90){
	    m_clsingletower->Fill((double)condarrcallo-condarrtkr);
	    m_clbytower[tkrvector[0]]->Fill((double)condarrcallo-condarrtkr);
	  }
	  if(calhivector.size()==1 && tkrvector[0]==calhivector[0]&& maxIdxface.getTwr()==calhivector[0] && maxeface>900 ){
	    m_chsingletower->Fill((double)condarrcalhi-condarrtkr);
	    m_chbytower[tkrvector[0]]->Fill((double)condarrcalhi-condarrtkr);
	  }
	}
	// only use event if max energy crystals is in the path of the track
	if (rec){
	  // extract track information from recon file
	  TkrRecon* tkrRecon = rec->getTkrRecon();
	  TObjArray* tracks = tkrRecon->getTrackCol();
	  if(tracks->GetEntries()>0){
	    //first track
	    TkrTrack* tkrTrack = dynamic_cast<TkrTrack*>(tracks->At(0));
	    if(tkrTrack) {
	      // end of track parameters
	      TkrTrackHit* hit = (TkrTrackHit*) tkrTrack->Last();
	      TVector3 endPos = hit->getPoint(TkrTrackHit::SMOOTHED);
	      TVector3 endSlopeTmp = hit->getDirection(TkrTrackHit::SMOOTHED);
	      TVector3 endSlope = endSlopeTmp.Unit();
	      // Need -1 here .....
	      const Vec3D tkr1EndDir(-1.0*endSlope.x(), -1.0*endSlope.y(),-1.0*endSlope.z());
	      double theta=tkr1EndDir.getTheta();
	      if (theta<MAXTHETARADS){
		const Vec3D tkr1EndPos(endPos.x(), endPos.y(), endPos.z());
		std::list<XtalIdx> tkrHitMap;
		tkrHitMap.clear();
		// find track intersection w/ each Cal lyr
		for (LyrNum lyr; lyr.isValid(); lyr++) {
		  // position of tracker track @ top of current lyr
		  const Vec3D   trkLyrTopPos(trkToZ (tkr1EndPos, tkr1EndDir,
						     // put me .001mm from z boundary to avoid
						     // gnarly floating point rounding issues
						     // most geometry is defined to .01 mm
						     lyrCtrZ (lyr)+CsIHeight/2-.001));
		  const XtalIdx xtalTop(pos2Xtal (trkLyrTopPos));

		  // make sure track passes through valid crystals
		  if (xtalTop.isValid()){
		    
		    const Vec3D   trkLyrCtrPos(trkToZ (tkr1EndPos, tkr1EndDir, lyrCtrZ (lyr)));
		    const XtalIdx xtalCtr(pos2Xtal (trkLyrCtrPos));

		    if (xtalCtr.isValid()){
		      const Vec3D   trkLyrBtmPos(trkToZ (tkr1EndPos, tkr1EndDir,
							 // put me .001 from z boundary to avoid
							 // gnarly floating point rounding issues
							 // most geometry is defined to .01 mm
							 lyrCtrZ (lyr)-CsIHeight/2+.001));
		      const XtalIdx xtalBtm(pos2Xtal (trkLyrBtmPos));

		      if (xtalBtm.isValid()){
			// track passes through valid crystals

			// check that entry / exit point are in same xtal
			// (guarantees entry & exit from same xtal top & bottom face
			if (xtalTop == xtalBtm){
			  // check that entry/exit point are in proper region of xtal face
			  const Vec3D ctrPos(xtalCtrPos (xtalCtr));

			  const DirNum dir(xtalCtr.getLyr ().getDir());
			  if (validXtalIntersect(dir, ctrPos, trkLyrTopPos) &&
			      validXtalIntersect(dir, ctrPos, trkLyrBtmPos)){

			    tkrHitMap.push_back(xtalCtr);
			  }
			}
		      }
		    }
		  }
		}
		
		// possible that there were no good xtal intersections
		if (!tkrHitMap.empty()){
		  if (maxeface>0 && std::find(tkrHitMap.begin(),tkrHitMap.end(),maxIdxface)!=tkrHitMap.end()){
		    if(gem.getConditionSummary()&0x4 && maxeface>90)m_clintersect->Fill((double)condarrcallo-condarrtkr);
		    if(gem.getConditionSummary()&0x8 && maxeface>900)m_chintersect->Fill((double)condarrcalhi-condarrtkr);
		  }
		  //		      std::cout<<counter<<" ";
		}
	      }
	    }
	  }
	}	

    }  // end analysis code in event loop
    

std::cout << "**** End of event loop **** " << std::endl;
std::cout << eng5<< " Engine 5 events" << std::endl;
std::cout << gamma<< " Gamma events" << std::endl;


}

void treqCAL::markError(char* line, TestReport *r){
  char text[128];
  m_status=1;
  strcpy(text,line);
  r->redtext(line,text); 
}
   
