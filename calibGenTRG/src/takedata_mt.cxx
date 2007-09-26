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
takedata_mt::takedata_mt() 
{
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

int takedata_mt::readParameterFile(char* pf){
  ifstream infile(pf);
  char line[8192];
  infile.getline(line,8000);
  map<string,string> params;
  map<string,string> properties;
  vector<string> keys;
  eval ev;
  map<string,string>::const_iterator it;
  params=ev.readdict(line);
  it=params.find("\'csvTestId\'");
  if (it==params.end())return -1;
  keys=ev.readlistortuple(it->second);
  m_testnr=ev.unquote(keys[1]);
  if (m_testnr!="pedestals")m_stepNumber=atoi(m_testnr.c_str());
  else return -1;
  cout<<"test id "<<m_stepNumber<<endl;
  for (int i=0;i<16;i++)m_temidlist[i]=i;
  it=params.find("\'runTimeOrig\'");
  if (it==params.end())return -1;
  m_evPerStep=(int)(atof(ev.unquote(it->second).c_str())*3600);
  cout<<"evPerStep"<<m_evPerStep<<endl;
  
  it=params.find("\'csvTestProperties\'");
  if (it==params.end())return -1;
  properties=ev.readdict(it->second);
  it=properties.find("(\'CAL\', \'TACK   delay [ticks]\')");
  if (it==properties.end())return -1;
  keys=ev.readlistortuple(it->second);
  m_tack_cal=atoi(ev.unquote(keys[0]).c_str());
  cout<<"CAL TACK "<<m_tack_cal<<endl;
  it=properties.find("(\'TKR\', \'TACK delay [ticks]\')");
  if (it==properties.end())return -1;
  keys=ev.readlistortuple(it->second);
  m_tack_tkr=atoi(ev.unquote(keys[0]).c_str());
  cout<<"TKR TACK "<<m_tack_tkr<<endl;
  it=properties.find("(\'ACD\', \'Hold delay [ticks]\')");
  if (it==properties.end())return -1;
  keys=ev.readlistortuple(it->second);
  m_tack_acd=atoi(ev.unquote(keys[0]).c_str());
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
    m_gid=evt->getMetaEvent().run().id();
    std::cout<<"gid "<<m_gid<<std::endl;
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
