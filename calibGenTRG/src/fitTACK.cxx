#include "fitTACK.h"
#include "TROOT.h"
#include "math.h"
#include "TestReport.h"

fittack::fittack(){
  first=0;
  last=0;
  npoints=0;
  useratiofile=false;
  
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
void fittack::setratiofilename(const char* name){
  useratiofile=true;
  for (int i=first;i<=last;i++){
    char fn[200];
    sprintf(fn,"%s_%d.root",name,i);
    ratiofilename[i]=fn;
  }
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
      if(!useratiofile){
	sprintf(histname,"calratio_tower_%d_step_%d",j,i);
	calratio[i][j]=(TH1D*)gDirectory->Get(histname);
	calratio[i][j]->SetDirectory(gROOT);
	if(j==0)calratiosum[i]=(TH1D*)calratio[i][j]->Clone();
	else calratiosum[i]->Add(calratio[i][j]);
      }
    }
    sprintf(histname,"adchist_alltowers_step_%d",i);
    sprintf(histtitle,"CAL ADC histogram all towers step %d",i);
    adchistsum[i]->SetNameTitle(histname,histtitle);
    sprintf(histname,"ACD_adchist_step_%d",i);
    acdadchist[i]=(TH1D*)gDirectory->Get(histname);
    acdadchist[i]->SetDirectory(gROOT);
    a.Close();
    if(useratiofile){
      a.Open(ratiofilename[i].c_str());
      for (int j=0;j<16;j++){
	sprintf(histname,"calratio_tower_%d_step_%d",j,i);
	calratio[i][j]=(TH1D*)gDirectory->Get(histname);
	calratio[i][j]->SetDirectory(gROOT);
	if(j==0)calratiosum[i]=(TH1D*)calratio[i][j]->Clone();
	else calratiosum[i]->Add(calratio[i][j]);
      }
      a.Close();
    }
    sprintf(histname,"calratio_alltowers_step_%d",i);
    sprintf(histtitle,"CAL HEX/LEX histogram all towers step %d",i);
    calratiosum[i]->SetNameTitle(histname,histtitle);
  }
}

 
void fittack::fit(){
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(11);
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
  TF1* eg=new TF1("eg",&egfun,35,70,3);
  //TF1* fnp=new TF1("fif","[3]+[4]*x+[5]*x*x+landau(0)" , 0,3);
  acdwaveform=new TH1D("acdwaveform","ACD waveform",npoints,tack_acd[first]-stepacd/2.,tack_acd[last]+stepacd/2.);
  langau->SetParameters(.1,1,30000,.15,1000);
  for (int i=first;i<=last;i++){
    fif->SetParameters(200000.,1.2,.2,10000,0.5);
    acdadchist[i]->Fit(langau,"em","",0.1,3);
    acdadchist[i]->Write();
    double meanfit=langau->GetParameter(1);
    double emeanfit=langau->GetParError(1);
    acdwaveform->SetBinContent(i-first+1,meanfit);
    acdwaveform->SetBinError(i-first+1,emeanfit);
  } 
  fn->SetParameters(1, 25,50,.5) ;
  acdwaveform->Fit(fn,"qem");
  acdwaveform->Write();
  acdwaveform->Draw();
  sprintf(epsname,"acd_%s.eps",textname.c_str());
  sprintf(gifname,"acd_%s.gif",textname.c_str());
  c1->SaveAs(epsname);
  c1->SaveAs(gifname);
  acdfitmean=fn->GetParameter(1);
  eacdfitmean=fn->GetParError(1);
  char histname[128],histtitle[128];
  calsumwaveform=new TH1D("calsumwaveform","CALLO Waveform All Towers",npoints,tack_cal[first]-stepcal/2.,tack_cal[last]+stepcal/2.);
  calhisumwaveform=new TH1D("calhisumwaveform","CALHI Waveform All Towers",npoints,tack_cal[first]-stepcal/2.,tack_cal[last]+stepcal/2.);
  calratiosumwaveform=new TH1D("calratiosumwaveform","CAL HEX/LEX Waveform All Towers",npoints,tack_cal[first]-stepcal/2.,tack_cal[last]+stepcal/2.);
  langau->SetParameters(0.06,11.3,10000,0.05,20);
  for (int i=0;i<16;i++){
    sprintf(histname,"calwaveform_tower_%d",i);
    sprintf(histtitle,"Cal Waveform Tower %d",i);
    calwaveform[i]=new TH1D(histname,histtitle,npoints,tack_cal[first]-stepcal/2.,tack_cal[last]+stepcal/2.);
    sprintf(histname,"calhiwaveform_tower_%d",i);
    sprintf(histtitle,"Calhi Waveform Tower %d",i);
    calhiwaveform[i]=new TH1D(histname,histtitle,npoints,tack_cal[first]-stepcal/2.,tack_cal[last]+stepcal/2.);
    sprintf(histname,"calratiowaveform_tower_%d",i);
    sprintf(histtitle,"Cal Ratio Waveform Tower %d",i);
    calratiowaveform[i]=new TH1D(histname,histtitle,npoints,tack_cal[first]-stepcal/2.,tack_cal[last]+stepcal/2.);
    for (int j=first;j<=last;j++){
      if (i==0){
	//	fifl->SetParameters(100000,15,2,0.5,5000,.1);
	//fnp->SetParameters(100000,15,2,5000,-100,0);
	//adchistsum[j]->Fit(fifl,"q","",10,30);
	adchistsum[j]->Fit(langau,"em","");
	//	adchistsum[j]->Fit(langau,"q","");
	adchistsum[j]->Write();
	calratiosum[j]->Write();
	double meanfit=langau->GetParameter(1);
	double emeanfit=langau->GetParError(1);
	calsumwaveform->SetBinContent(j-first+1,meanfit);
	calsumwaveform->SetBinError(j-first+1,emeanfit);
	double ratio=calratiosum[j]->GetMean();
	double nentries=calratiosum[j]->GetEntries();
	double eratio;
	if (nentries>0)
	  eratio=calratiosum[j]->GetRMS()/sqrt(float(nentries));
	else{
	  eratio=0;
	  ratio=0;
	}
	calratiosumwaveform->SetBinContent(j-first+1,ratio);
	calratiosumwaveform->SetBinError(j-first+1,eratio);
	calhisumwaveform->SetBinContent(j-first+1,meanfit*ratio);
	calhisumwaveform->SetBinError(j-first+1,sqrt(ratio*ratio*emeanfit*emeanfit+meanfit*meanfit*eratio*eratio));
      }
      fifl->SetParameters(7000,15,2,0.5,5000,.1);
      //fnp->SetParameters(100000,15,2,5000,-100,0);
      //adchist[j][i]->Fit(fifl,"q","",10,30);
      adchist[j][i]->Fit(langau,"em","");
      //adchist[j][i]->Fit(langau,"q","");
      adchist[j][i]->Write();
      calratio[j][i]->Write();
      double meanfit=langau->GetParameter(1);
      double emeanfit=langau->GetParError(1);
      calwaveform[i]->SetBinContent(j-first+1,meanfit);
      calwaveform[i]->SetBinError(j-first+1,emeanfit);
      double ratio=calratio[j][i]->GetMean();
      double nentries=calratio[j][i]->GetEntries();
      double eratio;
      if (nentries>0)
	eratio=calratio[j][i]->GetRMS()/sqrt(float(nentries));
      else{
	eratio=0;
	ratio=0;
      }
      calratiowaveform[i]->SetBinContent(j-first+1,ratio);
      calratiowaveform[i]->SetBinError(j-first+1,eratio);
      calhiwaveform[i]->SetBinContent(j-first+1,meanfit*ratio);
      calhiwaveform[i]->SetBinError(j-first+1,sqrt(ratio*ratio*emeanfit*emeanfit+meanfit*meanfit*eratio*eratio));

    }
    if (i==0){
      fn->SetParameters(10,50,80,.5) ;
      calsumwaveform->Fit(fn,"qem");
      //  eg->SetParameters(10,45,70);
       // calsumwaveform->Fit(eg,"q");
	calsumwaveform->Write();
	calsumwaveform->Draw();
	calfitmeansum=fn->GetParameter(1);
	ecalfitmeansum=fn->GetParError(1);
	sprintf(epsname,"cal_%s_sum.eps",textname.c_str());
	sprintf(gifname,"cal_%s_sum.gif",textname.c_str());
	c1->SaveAs(epsname);
	c1->SaveAs(gifname);
        eg->SetParameters(10,45,70);
      fn->SetParameters(10,50,80,.5) ;
      calhisumwaveform->Fit(fn,"qem");
      //  calhisumwaveform->Fit(eg,"q");
	calhisumwaveform->Write();
	calhisumwaveform->Draw();
	calhifitmeansum=fn->GetParameter(1);
	ecalhifitmeansum=fn->GetParError(1);
	sprintf(epsname,"calhi_%s_sum.eps",textname.c_str());
	sprintf(gifname,"calhi_%s_sum.gif",textname.c_str());
	c1->SaveAs(epsname);
	c1->SaveAs(gifname);
	calratiosumwaveform->Write();
    }
    fn->SetParameters(10,50,80,.5) ;
    calwaveform[i]->Fit(fn,"qem");
    //eg->SetParameters(10,45,70);
    //calwaveform[i]->Fit(eg,"q");
    calwaveform[i]->Write();
    calwaveform[i]->Draw();
    sprintf(epsname,"cal_%s_%d.eps",textname.c_str(),i);
    sprintf(gifname,"cal_%s_%d.gif",textname.c_str(),i);
    c1->SaveAs(epsname);
    c1->SaveAs(gifname);
    calfitmean[i]=fn->GetParameter(1);
    ecalfitmean[i]=fn->GetParError(1);
    //eg->SetParameters(10,45,70);
    fn->SetParameters(10,50,80,.5) ;
    calhiwaveform[i]->Fit(fn,"qem");
    //calhiwaveform[i]->Fit(eg,"q");
    calhiwaveform[i]->Write();
    calhiwaveform[i]->Draw();
    sprintf(epsname,"calhi_%s_%d.eps",textname.c_str(),i);
    sprintf(gifname,"calhi_%s_%d.gif",textname.c_str(),i);
    c1->SaveAs(epsname);
    c1->SaveAs(gifname);
    calhifitmean[i]=fn->GetParameter(1);
    ecalhifitmean[i]=fn->GetParError(1);
    calratiowaveform[i]->Write();
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
  else r.addstatus("Warning");
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
    sprintf(tem,"CALLO %d",i);
    sprintf(gifname,"cal_%s_%d.gif",textname.c_str(),i);
    r.addlink("Unit",gifname,tem); 
    sprintf(tem,"CALHI %d",i);
    sprintf(gifname,"calhi_%s_%d.gif",textname.c_str(),i);
    r.addlink("Unit",gifname,tem); 
    sprintf(tem,"TKR %d",i);
    sprintf(gifname,"tkr_%s_%d.gif",textname.c_str(),i);
    r.addlink("Unit",gifname,tem); 
  }
  sprintf(gifname,"cal_%s_sum.gif",textname.c_str());
  r.addlink("Unit",gifname,"CALLO All Towers");
  sprintf(gifname,"calhi_%s_sum.gif",textname.c_str());
  r.addlink("Unit",gifname,"CALHI All Towers");
  sprintf(gifname,"tkr_%s_sum.gif",textname.c_str());
  r.addlink("Unit",gifname,"TKR All Towers");
  sprintf(gifname,"acd_%s.gif",textname.c_str());
  r.addlink("Unit",gifname,"ACD");
  char* restable[]={"Subsystem","Status","Fit Delay+-Error","Best Delay Setting"};
  for (int i=0;i<16;i++){
    sprintf(name,"Results for tower at TEM %d",i);
    r.newheadline(name);
    r.starttable(restable,4);
    sprintf(line[0],"Calorimeter Low Energy");
    if (calstatus[i])r.redtext(line[1],"Warning");
    else r.greentext(line[1],"Passed");
    sprintf(line[2],"%.2f +- %.2f",calfitmean[i],ecalfitmean[i]);
    sprintf(line[3],"%d",int(calfitmean[i]+.5));
    r.addtableline(line,4);
    sprintf(line[0],"Calorimeter High Energy");
    if (calstatus[i])r.redtext(line[1],"Warning");
    else r.greentext(line[1],"Passed");
    sprintf(line[2],"%.2f +- %.2f",calhifitmean[i],ecalhifitmean[i]);
    sprintf(line[3],"%d",int(calhifitmean[i]+.5));
    r.addtableline(line,4);
    sprintf(line[0],"Tracker");
    if (tkrstatus[i])r.redtext(line[1],"Warning");
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
  sprintf(line[0],"Calorimeter Low Energy");
  if (calstatussum)r.redtext(line[1],"Warning");
  else r.greentext(line[1],"Passed");
  sprintf(line[2],"%.2f +- %.2f",calfitmeansum,ecalfitmeansum);
  sprintf(line[3],"%d",int(calfitmeansum+.5));
  r.addtableline(line,4);
  sprintf(line[0],"Calorimeter High Energy");
  if (calstatussum)r.redtext(line[1],"Warning");
  else r.greentext(line[1],"Passed");
  sprintf(line[2],"%.2f +- %.2f",calhifitmeansum,ecalhifitmeansum);
  sprintf(line[3],"%d",int(calhifitmeansum+.5));
  r.addtableline(line,4);
  sprintf(line[0],"Tracker");
  if (tkrstatussum)r.redtext(line[1],"Warning");
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
  if (acdstatus)r.redtext(line[1],"Warning");
  else r.greentext(line[1],"Passed");
  sprintf(line[2],"%.2f +- %.2f",acdfitmean,eacdfitmean);
  sprintf(line[3],"%d",int(acdfitmean+.5));
  r.addtableline(line,4);
  r.endtable();

    
  r.writereport();
}
  
  
