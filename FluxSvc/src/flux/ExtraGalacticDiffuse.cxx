// ExtraGalacticDiffuse.cxx: implementation of the ExtraGalacticDiffuse class.

//#include "StdAfx.h"
#include <cmath>
#include <iostream>
#include "ExtraGalacticDiffuse.h"
#include "time.h"
#include<vector>
#include "FluxSvc/SpectrumFactory.h"

static SpectrumFactory<ExtraGalacticDiffuse> factory;


long double ExtraGalacticDiffuse::pofi(long double intensity){  //this function gives P(I).  see documentation.
    long double p;
    //printf("\nabout to calculate pofi...");
    if(intensity>=pow(10.0,-7)){
        p=2.49555*pow(10.0,-13)*pow(intensity,-2.5);//egret range value for pofi
    }else{
        p=2.49555*pow(10.0,-13)*pow(intensity,-2.5-0.05*(7.0+(log(intensity)/log(10.0))));	
    }	//glast range value for pofi
    //printf("\np=%12.10e, intensity=%12.10e",p,intensity);
    
    return p;
}

//long double ExtraGalacticDiffuse::pofi(long double intensity);


long double ExtraGalacticDiffuse::random(){
    long double r=rand()/32767.0;
    return r;
}


ExtraGalacticDiffuse::DELTAX ExtraGalacticDiffuse::gaussianspread(){  //gives uncertainty to locations of incoming photons.
    DELTAX spread;
    long double x,y,theta,phi;
    
    theta=log(1.0/random());
    phi=random()*2.0*3.141592653589793238;
    
    x=theta*sin(phi);
    y=theta*cos(phi);
    
    spread.x=x;
    spread.y=y;
    return spread;
}


void ExtraGalacticDiffuse::addtophotons(long double x,long double y){
    int n;
    DELTAX smear;
    
    for(n=0 ; list[n].x!='\0' ; n++){
    }
    smear=gaussianspread();
    list[n].x=x+smear.x;
    list[n+1].x='\0';
    list[n].y=y+smear.y;
    list[n+1].y='\0';
    //	printf("\nAdded Photon to list, x=%lf, y=%lf",x,y);
}


void ExtraGalacticDiffuse::addtoctlg(long double x,long double y,long double intsty){
    int n;
    for(n=0 ; ctlg[n].x!='\0' ; n++){
    }
    ctlg[n].x=x;
    ctlg[n+1].x='\0';
    ctlg[n].y=y;
    ctlg[n].intensity=intsty;
}

void ExtraGalacticDiffuse::findandaddnew(){  //finds a new source to add from the catalog, adds it, 
    //and sends a photon from it.
    long double x=random()*skysizex;
    long double y=random()*skysizey;
    long double prob=random();
    //	cout << "x,y="<<  x  << '\n' << y  << '\n';
    long double dx=0.0000001;
    long double i=lowthresholdofintensity;
    
    while(prob > 0 && i<highthresholdofintensity){
        dx=0.01*pow(10.0,log10(i));
        prob-=(dx)*pofi(i);
        i+=dx;
        //	printf("\nin the findandaddnew loop; i=%12.10e, prob=%12.10e, dx=%12.10e , logi=%lf\n",i,prob,dx,log10(i));
    }
    
    addtophotons(x,y);
    addtoctlg(x,y,i-dx);
    RemInt-=i;
}

void ExtraGalacticDiffuse::sendphotonfromcatalog(){
    int i=0;
    long double num;
    long double catprob=1.0-NewProb;
    
    num=random()*catprob;
    for(i=0 ; num>0 ; i++){
        num-=ctlg[i].intensity/(TotInt-RemInt);}
    addtophotons(ctlg[i-1].x,ctlg[i-1].y);
}


ExtraGalacticDiffuse::PHOTON *ExtraGalacticDiffuse::create()
{
    long double rand;
    PHOTON *pointpho=NULL;
    
    int g=0;
    
    srand((unsigned)time(NULL));
    //					addtoctlg(0.5,2.9704,0.0001);
    //					RemInt-=0.0001;
    while(RemInt>=threshold || g<maxnumofphotons){
        //this sends a photon from the catalog if the random number was below the 
        //percentage of intensity already accounted for..
        g++;
        rand=random();
        NewProb=RemInt/TotInt;
        //	cout << "newprob=" << NewProb << '\n';
        if(rand<=NewProb){
            findandaddnew();
            //		printf("\ncompleted printandaddnew, Remint=%12.10e",RemInt);
        }
        else{sendphotonfromcatalog();  //and finds a new source if not.
        //printf("\ncompleted sendphotonfromcatalog");
        }
    }
    
    
    return pointpho;
};



//here's the constructor - grabs the data and makes a constant 

ExtraGalacticDiffuse::ExtraGalacticDiffuse(const char* name,float Emin, float Emax, float index){
    TotInt=.0000976759*multiplierduetosizeofsky;
    RemInt=TotInt;
    list[0].x='\0';
    list[0].y='\0';
    ctlg[0].x='\0';
    ctlg[0].y='\0';
    
    //Initialize the subclass with the appropriate values (from this object)
    SimpleSpectrum(name, Emin, Emax, index);
    
    //Call the create() function to initialize this object
    const ExtraGalacticDiffuse::PHOTON* photonlist=ExtraGalacticDiffuse::create();
    
    
    std::vector<std::pair<double,double> > listofphotons;
    
    
    //then put the photons in list into listofphotons
    
    int lenlst;
    
    std::pair<double,double> holder;
    //	cout << "totint=" << gms.TotInt << '\n';
    
    for(lenlst=0 ; list[lenlst].x!='\0' ; lenlst++){
        holder.first=list[lenlst].x;
        holder.second=list[lenlst].y;
        listofphotons.push_back (holder);
    }
    
    
    srcpnt = listofphotons.begin();
    //cout << "totint=" << TotInt << '\n';
    
};


//and our destructor:
ExtraGalacticDiffuse::~ExtraGalacticDiffuse(){
};


std::pair<double,double> ExtraGalacticDiffuse::dir(double e){
    
    std::pair<double,double> one;
    std::pair<double,double> coords;
    
    one=*srcpnt;
    
    coords.first=one.first;
    coords.second=one.second;
    srcpnt++;
    return coords;
}



