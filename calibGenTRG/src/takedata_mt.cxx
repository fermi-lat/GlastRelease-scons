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

#include "calibGenTRG/takedata_mt.h"
#include "eval.h"
#include "tackcuts.h"

#include "calibGenTRG/RootTreeAnalysis.h"
#include "TProfile.h"
#include "CalGeom.h"

#include <list>

using namespace CLHEP;
using namespace CalUtil;

ClassImp(takedata_mt)
  
// 
// Initialization
//
takedata_mt::takedata_mt() : m_caltuplefile(0),m_caltuple(0),m_merittuplefile(0), m_merittuple(0),m_useKalman(false), m_useToT(false),
			     m_kalmanLower(700), m_ToTLower(1), m_ToTUpper(1.6),
			     m_useOneTrack(false), m_useMip(false), m_MipLower(3),m_useptmaglat(0),m_ptmaglatlower(15),m_ptmaglatupper(25)
{
}

void takedata_mt::initCalTuple(const char* filename){
  if (m_caltuplefile){
    delete m_caltuple;
    delete m_caltuplefile;
  }
    m_caltuplefile=TFile::Open(filename, "READ");
    if(m_caltuplefile->IsOpen()){
      std::cout<<"        CAL:    "<<filename<<std::endl;
      m_caltuple=(TTree*)gDirectory->Get("CalTuple");
      m_caltuple->SetBranchAddress("CalXtalFaceSignalAllRange[16][8][12][2][4]",m_calrange);
    }else{
      std::cout<<"CAL file "<<filename<<" not found"<<std::endl;
    }
}
void takedata_mt::initMeritTuple(const char* filename){
  if (m_merittuplefile){
    delete m_merittuple;
    delete m_merittuplefile;
  }
    m_merittuplefile=TFile::Open(filename, "READ");
    if(m_merittuplefile->IsOpen()){
      std::cout<<"        Merit:    "<<filename<<std::endl;
      m_merittuple=(TTree*)gDirectory->Get("MeritTuple");
      m_merittuple->SetBranchStatus("*",0);
      m_merittuple->SetBranchStatus("PtMcIlwainL",1);
      m_merittuple->SetBranchAddress("PtMcIlwainL",&m_ptmaglat);
    }else{
      std::cout<<"Merit file "<<filename<<" not found"<<std::endl;
    }
}

    
void takedata_mt::inithistos(){

  char histname[128];
  char histtitle[128];
  for (int i=0;i<16;i++){
    m_totalcount[i]=m_totalsquare[i]=0;
    m_totalgood[i]=0;
  }
  for (int i=0;i<16;i++){
    sprintf(histname,"efficiency_tower_%d_step_%d" ,i,m_stepNumber);
    sprintf(histtitle,"Efficiency tower %d step %d", i,m_stepNumber);
    m_tkrocc[i]=new TH1D(histname,histtitle,101,-0.005,1.005);
    m_tkrocc[i]->GetXaxis()->SetTitle("Hit efficiency");

    sprintf(histname,"adchist_tower_%d_step_%d" ,i,m_stepNumber);
    sprintf(histtitle,"CAL ADC histogram tower %d step %d" ,i,m_stepNumber);
    m_adchist[i]= new TH1D(histname,histtitle,100,0,30);
    m_adchist[i]->SetMinimum(0.);
    m_adchist[i]->GetXaxis()->SetTitle("MeV");
    sprintf(histname,"calratio_tower_%d_step_%d" ,i,m_stepNumber);
    sprintf(histtitle,"HEX/LEX ratio in tower %d step %d" ,i,m_stepNumber);
    m_calratio[i]= new TH1D(histname,histtitle,151,0.495,2.005);
    m_calratio[i]->SetMinimum(0.);
    m_calratio[i]->GetXaxis()->SetTitle("Ratio");

  }
  sprintf(histname,"ACD_adchist_step_%d" ,m_stepNumber);
  sprintf(histtitle,"ACD ADC histogram step %d" ,m_stepNumber);
  m_acdadchist=new TH1D(histname,histtitle,100,0,3);
  m_acdadchist->GetXaxis()->SetTitle("MIPS");
}  

void takedata_mt::savehistos(char* filename){
  strcpy(m_filename,filename);
  TFile f(filename,"recreate");
  for (int i=0;i<16;i++){
    m_tkrocc[i]->Write();
    m_adchist[i]->Write();
    m_calratio[i]->Write();
  }
  m_acdadchist->Write();
  f.Close();
}
  
void takedata_mt::writeoutresults(char* rf){
  ofstream res(rf);
  res<<m_gid<<std::endl;
  res<<m_evPerStep<<std::endl;
  res<<m_nev<<std::endl;
  res<<m_tack_cal<<std::endl;
  res<<m_tack_acd<<std::endl;
  res<<m_tack_tkr<<std::endl;
  res<<m_filename<<std::endl;
  for (int i=0;i<16;i++){
    res<<m_totalcount[i]<<std::endl;
    res<<m_totalgood[i]<<std::endl;
    res<<m_totalsquare[i]<<std::endl;
  }
  
}

int takedata_mt::setParameters(int step, int caldelay, int tkrdelay, int acddelay){
  m_stepNumber=step;
  cout<<"Step number "<<m_stepNumber<<std::endl;
  m_tack_cal=caldelay;
  cout<<"CAL TACK "<<m_tack_cal<<endl;
  m_tack_tkr=tkrdelay;
  cout<<"TKR TACK "<<m_tack_tkr<<endl;
  m_tack_acd=acddelay;
  cout<<"ACD TACK "<<m_tack_acd<<endl;
  return 0;
}


//
// Main event loop
//
void takedata_mt::Go(Long64_t numEvents)
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
    m_gid=evt->getRunId();
    std::cout<<"gid "<<m_gid<<std::endl;
    m_nev=nEvtMax-m_StartEvent;
    m_evPerStep=m_nev;
    cout<<"Events used: "<<m_evPerStep<<endl;
    //}
    // BEGINNING OF EVENT LOOP
    for (Long64_t ievent=m_StartEvent; ievent<nEvtMax; ievent++) {
        if(ievent%10000==0) 
            std::cout << "\r** Processing Event " << ievent << std::flush;  
        
        if (mc) mc->Clear();
        if (evt) evt->Clear();
        if (rec) rec->Clear();
	if (m_merittuple)m_merittuple->GetEvent(ievent);
	if (m_useptmaglat&& (m_ptmaglat<m_ptmaglatlower || m_ptmaglat>m_ptmaglatupper))continue;
        Int_t nb = GetEvent(ievent);
	if (m_caltuple)m_caltuple->GetEvent(ievent);
	const Gem gem=evt->getGem();       

	int miplayers[8];
	double Elayers[8];
	for (int i=0;i<8;i++){
	  miplayers[i]=0;
	  Elayers[i]  =0;
	  }
	XtalIdx maxIdxface;

	if (rec){
	  if (evt->getL1T().getGemEngine()==4){ // Why ?
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
		    
		    // if the cristal energy is high enough, calculate the ratio
		    if (xenergy>100&&xenergy<800){
		      if (m_calrange[tower][layer][column][0][1]>0 && m_calrange[tower][layer][column][0][2]>0
			  && m_calrange[tower][layer][column][1][1]>0 && m_calrange[tower][layer][column][1][2]>0){
			float faceratio2=m_calrange[tower][layer][column][0][2]/m_calrange[tower][layer][column][1][2];
			float faceratio1=m_calrange[tower][layer][column][0][1]/m_calrange[tower][layer][column][1][1];
			if(faceratio1<1)faceratio1=1/faceratio1;
			if(faceratio2<1)faceratio2=1/faceratio2;
			if (faceratio1<2 && faceratio2<2){
			  float ratio1=m_calrange[tower][layer][column][0][2]/m_calrange[tower][layer][column][0][1];
			  float ratio2=m_calrange[tower][layer][column][1][2]/m_calrange[tower][layer][column][1][1];
			  m_calratio[tower]->Fill(ratio1);
			  m_calratio[tower]->Fill(ratio2);
			}
		      }
		    }
		  } //end of while	      
		}
	      }
	    } // end of getGemEngine==4

	    // Add here MIP selection from the calorimeter     
	  if (m_useMip){
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
	  	  Elayers[layer]+=xenergy; // and calculate Elayer for the whole LAT
		  miplayers[layer]++;
	  	  
	  	} //end of while	
	      }
	    }
	  } // end of if use mip

	    
	  // Check Layers Energy to see whether there is a mip in there or not
	  // layer counts for all CAL towers, this in someway cleans more events than checking only one tower.
	  int nlayers=0;
	  for (int i=0;i<8;i++){
	    if (Elayers[i]>5&&Elayers[i]<25) nlayers++;  // Require a mip in a single xtal
	  }
	  bool occveto=false; 
	  for (int i=0;i<8;i++){
	    if (miplayers[i]>2)occveto=true;
	  }
	  if(occveto)continue;
	  // extract track information from recon file
	  TkrRecon* tkrRecon = rec->getTkrRecon();
	  TObjArray* tracks = tkrRecon->getTrackCol();
	  if (tracks->GetEntries()==0)continue;
	  if(tracks->GetEntries()>0){

  	    // require TKR and ROI exactly
	    if( !(gem.getConditionSummary()&0x3) ) continue;  

	    // Require a MIP track in the calorimeter
	    if(m_useMip && nlayers<m_MipLower)continue;

	    // Select here event that have more than 1 track
	    if( (tracks->GetEntries())!=1 && m_useOneTrack)continue;

	    //first track
	    TkrTrack* tkrTrack = dynamic_cast<TkrTrack*>(tracks->At(0));
	    if(tkrTrack) {
	      // event cut on Kalman energy 
	      if(m_useKalman && tkrTrack->getKalEnergy( )<m_kalmanLower)continue;
	      // event cut track chi square 
	      if( (tkrTrack->getChiSquareSmooth()>20) && m_useOneTrack)continue; 


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
		  }
		  
		  Tkr_1_ToTAve += mips;
		}
		
		if(clustersOnTrack>3) {
		  Tkr_1_ToTAve /= clustersOnTrack;
		}
		if (Tkr_1_ToTAve>m_ToTUpper || Tkr_1_ToTAve<m_ToTLower)continue; // only use ToT between 1 and 1.6
	      }
		
	      
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
		  TObjArray *xtalRecCol = rec->getCalRecon()->getCalXtalRecCol();
		  TIter xtalIter(xtalRecCol);
		  CalXtalRecData *xtal = 0;
		  if (xtalRecCol){
		    if(xtalRecCol->GetEntries()>0 && xtalRecCol->GetEntries()<MAXCALHITS){
		      unsigned counter=0;
		      while ((xtal = (CalXtalRecData*)xtalIter.Next())) {
			const idents::CalXtalId cId(xtal->getPackedId());
			Int_t tower  = cId.getTower();
			const XtalIdx xtalIdx(cId);
			if (std::find(tkrHitMap.begin(),tkrHitMap.end(),xtalIdx)!=tkrHitMap.end()){
			  double xenergy=xtal->getEnergy()*cos(theta);
			  m_adchist[tower]->Fill(xenergy);
			  counter++;
			}
		      }
		      //		      std::cout<<counter<<" ";
		    }
		  }
		}
	      }
	    }
	  }

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
	      }
	    }	
	  }  
	}
	double eff[16];
	int nlayers[16];
	int tkrhits[16][36];
	for (int i=0;i<16;i++){
	  eff[i]=0;
	  nlayers[i]=0;
	  for (int j=0;j<36;j++){
	    tkrhits[i][j]=0;
	  }
	}
	const TObjArray* tkrDigiCol = evt->getTkrDigiCol();
	if (tkrDigiCol){
	  TIter tkrIter(tkrDigiCol);
	  TkrDigi *tkr = 0;
	  while ((tkr = (TkrDigi*)tkrIter.Next())) {
	    Int_t tower   = tkr->getTower().id();
	    Int_t bilayer = tkr->getBilayer();
	    Int_t view = tkr->getView();
	    if (bilayer%2) tkrhits[tower][bilayer*2+view]=1;
	    else tkrhits[tower][bilayer*2+1-view]=1;
	  }
	  for (int to=0;to<16;to++){
	    int first=35;
	    int last=0;
	    for (int i=0;i<36;i+=2){
	      if(tkrhits[to][i] && tkrhits[to][i+1])last=i;
	    }
	    for (int i=34;i>-1;i-=2){
	      if (tkrhits[to][i] && tkrhits[to][i+1])first=i;
	    }
	    if (first<last){
	      int count=0;
	      int goodcount=0;
	      for (int i=first+2;i<last;i++){
		count++;
		if (tkrhits[to][i]){
		  goodcount++;
		  m_totalsquare[to]++;
		}
	      }
	      if (count>0){
		eff[to]=double(goodcount)/double(count);
		nlayers[to]=count;
		m_totalgood[to]+=eff[to];
		m_totalcount[to]++;
		m_tkrocc[to]->Fill(eff[to],nlayers[to]);
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
