// ExtraGalacticDiffuse.h: interface for the ExtraGalacticDiffuse class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_GAMMASOURCE_H__A03408DA_6202_4ADF_8F8E_AFFF184A7AC2__INCLUDED_)
#define AFX_GAMMASOURCE_H__A03408DA_6202_4ADF_8F8E_AFFF184A7AC2__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

//#pragma inline_depth (2)
//Code for producing photon sets
#include "FluxSvc/SimpleSpectrum.h"
#include <cmath>
#include <algorithm>
#include <functional>
#include "CLHEP/Random/Random.h"
#include "FluxSvc/GPS.h"
#include <vector>

/*! Currently unused class simulating diffuse particle background.
*/
class ExtraGalacticDiffuse : public SimpleSpectrum {
public:
    //      protected:
    
#define maxnumofphotons 500  //maximum number of photons to ever send (the program continues sending even after all intensity is accounted for.
#define lowthresholdofintensity 0.000000001  //defines the minimum source intensity considered
#define highthresholdofintensity .0001 //for debugging and checking purposes
#define threshold 0.000000001  //for determining when all intensity is accounted for.
#define skysizex 10.0
#define skysizey 10.0
#define multiplierduetosizeofsky 1.0/(36.0*36.0)
    ///TotInt=.0000976759*multiplierduetosizeofsky;  
    ///you might want to change this to the actual calculated value 
    ///- it is in photons/sec/cm^2 over the whole sky
    long double TotInt;
    
    long double RemInt;
    //RemInt=TotInt;
    long double NewProb;
    
    
    typedef struct{
        long double x;
        long double y;
    }PHOTON;  //photon coodinates
    
    PHOTON list[50000];
    
    typedef struct{
        long double x;
        long double y;
        long double intensity;
    }SOURCE;  //format for a luminous source
    
    SOURCE ctlg[5000];
    
    
    typedef struct{
        long double x;
        long double y;
    }DELTAX;  //for adding gaussian smear to data
    
    
    std::vector<std::pair<double,double> >::iterator srcpnt;
    
    long double pofi(long double intensity);
    
    long double random();
    
    DELTAX gaussianspread();
    
    void addtophotons(long double x,long double y);
    
    void addtoctlg(long double x,long double y,long double intsty);
    
    void findandaddnew();
    
    void sendphotonfromcatalog();
    
    PHOTON *create();
    
    ~ExtraGalacticDiffuse();
    
    ///default constructor
    ExtraGalacticDiffuse(){}
    
    ExtraGalacticDiffuse(const char* name,float Emin, float Emax, float index);
    
    std::pair<double,double> dir(double e);
    
    
}; //class ending bracket

#endif // !defined(AFX_GAMMASOURCE_H__A03408DA_6202_4ADF_8F8E_AFFF184A7AC2__INCLUDED_)
