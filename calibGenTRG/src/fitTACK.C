#include <TH1D.h>
#include <TF1.h>
#include <fstream>
#include <iostream>
#include "TestReport.h"
#include <TStyle.h>
#include <TCanvas.h>
#include <time.h>
#include "LangauFun.h"

#define MAXCAL 5
#define MAXTKR 0
#define MAXACD 5
#define STDACD 24
#define STDCAL 45
#define STDTKR 0
using namespace std;

class fittack{
public:
  fittack();
  void settextinput(char*);
  void setrootoutput(char*);
  void setfirstlastpoint(int,int);
  void readtextfiles();
  void readrootfiles();
  TH1D* calhist(int i, int j){return adchist[i][j];}
  TH1D* acdhist(int i){return acdadchist[i];}
  TH1D* getAcdWaveform(){return acdwaveform;}
  void fit();
  int validate();
  void writereport(int status);
private:
  string textname;
  string rootoutput;
  int first,last,npoints;
  int evperstep[30],nev[30],tack_cal[30],tack_tkr[30],tack_acd[30];
  char gid[30][128];
  double stepacd,stepcal,steptkr;
  int totalcount[30][16],totalsquare[30][16];
  double totalgood[30][16];
  int totalcountsum[30],totalsquaresum[30];
  double totalgoodsum[30];
  string filename[30];
  double calfitmeansum,ecalfitmeansum;
  double tkrfitmeansum,etkrfitmeansum;
  double calfitmean[16],ecalfitmean[16];
  double tkrfitmean[16],etkrfitmean[16];
  double acdfitmean,eacdfitmean;
  int acdstatus,calstatus[16],tkrstatus[16];
  int calstatussum,tkrstatussum;
  TH1D *tkrocc[30][16], *adchist[30][16], *acdadchist[30];
  TH1D  *adchistsum[30];
  TH1D *acdwaveform, *calwaveform[16], *tkrwaveform[30];
  TH1D *calsumwaveform,*tkrsumwaveform;
};

fittack::fittack(){
  first=0;
  last=0;
  npoints=0;
  
}
void fittack::settextinput(char* name){
  textname=name;
}
void fittack::setrootoutput(char* name){
  rootoutput=name;
}
void fittack::setfirstlastpoint(int i,int j){
  first=i;
  last=j;
  npoints=last-first+1;
}
void fittack::readtextfiles(){
  for (int i=first;i<=last;i++){
    char name[128];
    sprintf(name,"%s_%d.txt",textname.c_str(),i);
    std::ifstream ifs;
    ifs.open(name);
    ifs>>gid[i]>>evperstep[i]>>nev[i]>>tack_cal[i]>>tack_acd[i]>>tack_tkr[i]>>filename[i];
    for (int j=0;j<16;j++){
      ifs>>totalcount[i][j];
      ifs>>totalgood[i][j];
      ifs>>totalsquare[i][j];
    }
  }
  if (first!=last){
    stepacd=(double)(tack_acd[first+1]-tack_acd[first]);
    stepcal=(double)(tack_cal[first+1]-tack_cal[first]);
    steptkr=(double)(tack_tkr[first+1]-tack_tkr[first]);
  }else{
    stepacd=5;
    steptkr=5;
    stepcal=5;
  }
}
  
void fittack::readrootfiles(){
  TFile a;
  char histname[127],histtitle[127];
  for (int i=first;i<=last;i++){
    a.Open(filename[i].c_str());
    for (int j=0;j<16;j++){
      sprintf(histname,"efficiency_tower_%d_step_%d",j,i);
      tkrocc[i][j]=(TH1D*)gDirectory->Get(histname);
      tkrocc[i][j]->SetDirectory(gROOT);
      sprintf(histname,"adchist_tower_%d_step_%d",j,i);
      adchist[i][j]=(TH1D*)gDirectory->Get(histname);
      adchist[i][j]->SetDirectory(gROOT);
      if(j==0)adchistsum[i]=(TH1D*)adchist[i][j]->Clone();
      else adchistsum[i]->Add(adchist[i][j]);
    }
    sprintf(histname,"adchist_alltowers_step_%d",i);
    sprintf(histtitle,"CAL ADC histogram all towers step %d",i);
    adchistsum[i]->SetNameTitle(histname,histtitle);
    sprintf(histname,"ACD_adchist_step_%d",i);
    acdadchist[i]=(TH1D*)gDirectory->Get(histname);
    acdadchist[i]->SetDirectory(gROOT);
    a.Close();
  }
}

 
void fittack::fit(){
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  gStyle->SetErrorX(0.002);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(.8);
  TCanvas *c1=new TCanvas();
  c1->Draw();
  char epsname[128];
  char gifname[128];
  TFile a(rootoutput.c_str(),"recreate");
  TF1* langau=new TF1(LangauFun::getLangauDAC());
  TF1* fifl=new TF1("fifl","landau" , 0,3);
  TF1* fif=new TF1("fif","[3]*exp(-x*[4])+landau(0)" , 0,3);
  TF1* fn=new TF1("fn","[0]*exp(-0.5*(pow(log(1+[3]*(x-[1])/[2]*sinh([3]*1.177)/1.177/[3])/[3],2)+[3]*[3]))",0,700);
  TF1* fnp=new TF1("fif","[3]+[4]*x+[5]*x*x+landau(0)" , 0,3);
  acdwaveform=new TH1D("acdwaveform","ACD waveform",npoints,tack_acd[first]-stepacd/2.,tack_acd[last]+stepacd/2.);
  for (int i=first;i<=last;i++){
    fif->SetParameters(200000.,1.2,.2,10000,0.5);
    acdadchist[i]->Fit(fif,"q","",0.5,3);
    acdadchist[i]->Write();
    double meanfit=fif->GetParameter(1);
    double emeanfit=fif->GetParError(1);
    acdwaveform->SetBinContent(i-first+1,meanfit);
    acdwaveform->SetBinError(i-first+1,emeanfit);
  } 
  fn->SetParameters(1.5, 25,5,.1) ;
  acdwaveform->Fit(fn,"q");
  acdwaveform->Write();
  acdwaveform->Draw();
  sprintf(epsname,"acd_%s.eps",textname.c_str());
  sprintf(gifname,"acd_%s.gif",textname.c_str());
  c1->SaveAs(epsname);
  c1->SaveAs(gifname);
  acdfitmean=fn->GetParameter(1);
  eacdfitmean=fn->GetParError(1);
  char histname[128],histtitle[128];
  calsumwaveform=new TH1D("calsumwaveform","CAL Waveform All Towers",npoints,tack_cal[first]-stepcal/2.,tack_cal[last]+stepcal/2.);
  for (int i=0;i<16;i++){
    sprintf(histname,"calwaveform_tower_%d",i);
    sprintf(histtitle,"Cal Waveform Tower %d",i);
    calwaveform[i]=new TH1D(histname,histtitle,npoints,tack_cal[first]-stepcal/2.,tack_cal[last]+stepcal/2.);
    for (int j=first;j<=last;j++){
      if (i==0){
	//	fifl->SetParameters(100000,15,2,0.5,5000,.1);
	//fnp->SetParameters(100000,15,2,5000,-100,0);
	//adchistsum[j]->Fit(fifl,"q","",10,30);
	adchistsum[j]->Fit(langau,"q","");
	adchistsum[j]->Write();
	double meanfit=fifl->GetParameter(1);
	double emeanfit=fifl->GetParError(1);
	calsumwaveform->SetBinContent(j-first+1,meanfit);
	calsumwaveform->SetBinError(j-first+1,emeanfit);
      }
      fifl->SetParameters(7000,15,2,0.5,5000,.1);
      //fnp->SetParameters(100000,15,2,5000,-100,0);
      //adchist[j][i]->Fit(fifl,"q","",10,30);
      adchist[j][i]->Fit(langau,"q","");
      adchist[j][i]->Write();
      double meanfit=fifl->GetParameter(1);
      double emeanfit=fifl->GetParError(1);
      calwaveform[i]->SetBinContent(j-first+1,meanfit);
      calwaveform[i]->SetBinError(j-first+1,emeanfit);
    }
    if (i==0){
        fn->SetParameters(12,45,100,.5) ;
	calsumwaveform->Fit(fn,"q");
	calsumwaveform->Fit(fn,"q");
	calsumwaveform->Fit(fn,"q");
	calsumwaveform->Write();
	calsumwaveform->Draw();
	calfitmeansum=fn->GetParameter(1);
	ecalfitmeansum=fn->GetParError(1);
	sprintf(epsname,"cal_%s_sum.eps",textname.c_str());
	sprintf(gifname,"cal_%s_sum.gif",textname.c_str());
	c1->SaveAs(epsname);
	c1->SaveAs(gifname);
    }
    fn->SetParameters(12,45,100,.5) ;
    calwaveform[i]->Fit(fn,"q");
    calwaveform[i]->Fit(fn,"q");
    calwaveform[i]->Fit(fn,"q");
    calwaveform[i]->Write();
    calwaveform[i]->Draw();
    sprintf(epsname,"cal_%s_%d.eps",textname.c_str(),i);
    sprintf(gifname,"cal_%s_%d.gif",textname.c_str(),i);
    c1->SaveAs(epsname);
    c1->SaveAs(gifname);
    calfitmean[i]=fn->GetParameter(1);
    ecalfitmean[i]=fn->GetParError(1);
  }
  tkrsumwaveform=new TH1D("tkrsumwaveform","TKR Waveform All Towers",npoints,tack_tkr[first]-steptkr/2.,tack_tkr[last]+steptkr/2.);
  for (int j=first;j<=last;j++){
    totalgoodsum[j]=0;
    totalcountsum[j]=0;
    totalsquaresum[j]=0;
  }
  for (int i=0;i<16;i++){
    sprintf(histname,"tkrwaveform_tower_%d",i);
    sprintf(histtitle,"Tkr Waveform Tower %d",i);
    tkrwaveform[i]=new TH1D(histname,histtitle,npoints,tack_tkr[first]-steptkr/2.,tack_tkr[last]+steptkr/2.);
    for (int j=first;j<=last;j++){
      tkrocc[j][i]->Write();
      totalgoodsum[j]+=totalgood[j][i];
      totalcountsum[j]+=totalcount[j][i];
      totalsquaresum[j]+=totalsquare[j][i];
      if(totalcount[j][i]>0){
	double tme2=totalgood[j][i]/double(totalcount[j][i]);
	double tmee2=sqrt(tme2*(1-tme2)/totalcount[j][i]);
	tkrwaveform[i]->SetBinContent(j-first+1,tme2);
	tkrwaveform[i]->SetBinError(j-first+1,tmee2);
      }
    }
    double mnval=10000;
    for (int j=first;j<=last;j++){
      double cn=tkrwaveform[i]->GetBinContent(j-first+1);
      if (cn!=0 && cn <mnval)mnval=cn;
    }
    tkrwaveform[i]->SetMinimum(mnval*.99);
    fn->SetParameters(10,0,20,1);
    tkrwaveform[i]->Fit(fn,"q");
    tkrwaveform[i]->Write();
    tkrwaveform[i]->Draw();
    sprintf(epsname,"tkr_%s_%d.eps",textname.c_str(),i);
    sprintf(gifname,"tkr_%s_%d.gif",textname.c_str(),i);
    c1->SaveAs(epsname);
    c1->SaveAs(gifname);
    tkrfitmean[i]=fn->GetParameter(1);
    etkrfitmean[i]=fn->GetParError(1);
  }
  for (int j=first;j<=last;j++){
    if(totalcountsum[j]>0){
      cout<<totalgoodsum[j]<<" "<<totalcountsum[j]<<endl;
      double tme2=totalgoodsum[j]/double(totalcountsum[j]);
      double tmee2=sqrt(tme2*(1-tme2)/totalcountsum[j]);
      tkrsumwaveform->SetBinContent(j-first+1,tme2);
      tkrsumwaveform->SetBinError(j-first+1,tmee2);
    } 
  }
  fn->SetParameters(10,0,20,1);
  tkrsumwaveform->Fit(fn,"0q");
  tkrfitmeansum=fn->GetParameter(1);
  etkrfitmeansum=fn->GetParError(1);
  double mnval=10000;
  for (int j=first;j<=last;j++){
    double cn=tkrsumwaveform->GetBinContent(j-first+1);
    if (cn!=0 && cn <mnval)mnval=cn;
  }
  tkrsumwaveform->SetMinimum(mnval*.99);
  tkrsumwaveform->Draw();
  tkrsumwaveform->Write();
  sprintf(epsname,"tkr_%s_sum.eps",textname.c_str());
  sprintf(gifname,"tkr_%s_sum.gif",textname.c_str());
  c1->SaveAs(epsname);
  c1->SaveAs(gifname);
  a.Close();
}

int fittack::validate(){
  int status=0;
  if (ecalfitmeansum>2 || fabs(calfitmeansum-STDCAL)>MAXCAL){
    status=1;
    calstatussum=1;
  }else  calstatussum=0;
  double bestdelay=tkrfitmeansum;
  if (bestdelay<0)bestdelay=0;
  if (fabs(bestdelay-STDTKR)>MAXTKR){
    tkrstatussum=1;
    status=1;
  }else tkrstatussum=0;
  for (int i=0;i<16;i++){
    if (ecalfitmean[i]>2 || fabs(calfitmean[i]-STDCAL)>MAXCAL){
      status=1;
      calstatus[i]=1;
    }else calstatus[i]=0;
    double bestdelay=tkrfitmean[i];
    if (bestdelay<0)bestdelay=0;
    if (fabs(bestdelay-STDTKR)>MAXTKR){
      tkrstatus[i]=1;
      status=1;
    }else tkrstatus[i]=0;
  }
  if (eacdfitmean>2 || fabs(acdfitmean-STDACD)>MAXACD){
    status=1;
    acdstatus=1;
  }else acdstatus=0;
  return status;
}

void fittack::writereport(int s){
  char name[128];
  sprintf(name,"report_%s.html",textname.c_str());
  TestReport r(name);
  r.newheadline("Test identification");
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  sprintf(name,"%s",asctime(timeinfo));
  r.additem("Date", name);
  r.additem("Application", "fitTACK.C");
  r.additem("Description", "TACK delay time-in analysis for multiple towers and ACD");
  r.additem("Type of run","Offline");
  if (s==0)r.addstatus("Passed");
  else r.addstatus("Failed");
  r.newheadline("Runs");
  char* header[]={"Ground Id","Number of events","Number of seconds","Trigger source","CAL TACK setting","TKR TACK setting","ACD TACK setting"};
  r.starttable(header,7);
  char* line[7];
  for (int j=0;j<7;j++)line[j]=new char[64];
  for (int i=first; i<=last;i++){
    sprintf(line[0],"%s",gid[i]);
    sprintf(line[1],"%d",nev[i]);
    sprintf(line[2],"%d",evperstep[i]);
    sprintf(line[3],"TKR");
    sprintf(line[4],"%d",tack_cal[i]);
    sprintf(line[5],"%d",tack_tkr[i]);
    sprintf(line[6],"%d",tack_acd[i]);
    r.addtableline(line,7);
  }
  r.endtable();
  r.newheadline("Plots");
  char tem[12],gifname[128];
  for (int i=0;i<16;i++){
    sprintf(tem,"CAL %d",i);
    sprintf(gifname,"cal_%s_%d.gif",textname.c_str(),i);
    r.addlink("Unit",gifname,tem); 
    sprintf(tem,"TKR %d",i);
    sprintf(gifname,"tkr_%s_%d.gif",textname.c_str(),i);
    r.addlink("Unit",gifname,tem); 
  }
  sprintf(gifname,"cal_%s_sum.gif",textname.c_str());
  r.addlink("Unit",gifname,"CAL All Towers");
  sprintf(gifname,"tkr_%s_sum.gif",textname.c_str());
  r.addlink("Unit",gifname,"TKR All Towers");
  sprintf(gifname,"acd_%s.gif",textname.c_str());
  r.addlink("Unit",gifname,"ACD");
  char* restable[]={"Subsystem","Status","Fit Delay+-Error","Best Delay Setting"};
  for (int i=0;i<16;i++){
    sprintf(name,"Results for tower at TEM %d",i);
    r.newheadline(name);
    r.starttable(restable,4);
    sprintf(line[0],"Calorimeter");
    if (calstatus[i])r.redtext(line[1],"Failed");
    else r.greentext(line[1],"Passed");
    sprintf(line[2],"%.2f +- %.2f",calfitmean[i],ecalfitmean[i]);
    sprintf(line[3],"%d",int(calfitmean[i]+.5));
    r.addtableline(line,4);
    sprintf(line[0],"Tracker");
    if (tkrstatus[i])r.redtext(line[1],"Failed");
    else r.greentext(line[1],"Passed");
    sprintf(line[2],"%.2f +- %.2f",tkrfitmean[i],etkrfitmean[i]);
    int bestvalue=int(tkrfitmean[i]+.5);
    if (bestvalue<0)bestvalue=0;
    sprintf(line[3],"%d",bestvalue);
    r.addtableline(line,4);
    r.endtable();
  }
  r.newheadline("Results for all 16 towers");
  r.starttable(restable,4);
  sprintf(line[0],"Calorimeter");
  if (calstatussum)r.redtext(line[1],"Failed");
  else r.greentext(line[1],"Passed");
  sprintf(line[2],"%.2f +- %.2f",calfitmeansum,ecalfitmeansum);
  sprintf(line[3],"%d",int(calfitmeansum+.5));
  r.addtableline(line,4);
  sprintf(line[0],"Tracker");
  if (tkrstatussum)r.redtext(line[1],"Failed");
  else r.greentext(line[1],"Passed");
  sprintf(line[2],"%.2f +- %.2f",tkrfitmeansum,etkrfitmeansum);
  int bestvalue=int(tkrfitmeansum+.5);
  if (bestvalue<0)bestvalue=0;
  sprintf(line[3],"%d",bestvalue);
  r.addtableline(line,4);
  r.endtable();
  r.newheadline("Results for ACD");
  r.starttable(restable,4);
  sprintf(line[0],"ACD");
  if (acdstatus)r.redtext(line[1],"Failed");
  else r.greentext(line[1],"Passed");
  sprintf(line[2],"%.2f +- %.2f",acdfitmean,eacdfitmean);
  sprintf(line[3],"%d",int(acdfitmean+.5));
  r.addtableline(line,4);
  r.endtable();

    
  r.writereport();
}
  
  
