//  The below code is a modified version of the ROOT/TMinuit class

// in order to have a stand alone C++ version of the Minuit package.

//

//  The usage of TString has been replaced by the simpler class 

// MString that has an STL API. The reference to TNamed, 

// TObject and TMethods has been removed. In ROOT it could be 

// reintroduced through a class TMinuit that could inherit the 

// Midnight one.

//

//    G.Barrand



/*************************************************************************

 * Copyright(c) 1995-1998, The ROOT System, All rights reserved.         *

 * Authors: Rene Brun, Nenad Buncic, Valery Fine, Fons Rademakers.       *

 *                                                                       *

 * Permission to use, copy, modify and distribute this software and its  *

 * documentation for non-commercial purposes is hereby granted without   *

 * fee, provided that the above copyright notice appears in all copies   *

 * and that both the copyright notice and this permission notice appear  *

 * in the supporting documentation. The authors make no claims about the *

 * suitability of this software for any purpose.                         *

 * It is provided "as is" without express or implied warranty.           *

 *************************************************************************/



//*-*-*-*-*-*-*-*-*-*-*-*The Minimization package*-*--*-*-*-*-*-*-*-*-*-*-*

//*-*                    ========================                         *

//*-*                                                                     *

//*-*   This package was originally written in Fortran by Fred James      *

//*-*   and part of PACKLIB (patch D506)                                  *

//*-*                                                                     *

//*-*   It has been converted to a C++ class  by R.Brun                   *

//*-*   The current implementation in C++ is a straightforward conversion *

//*-*   of the original Fortran version: The main changes are:            *

//*-*                                                                     *

//*-*   - The variables in the various Minuit labelled common blocks      *

//*-*     have been changed to the Minuit class data members.            *

//*-*   - The include file Minuit.h has been commented as much as possible*

//*-*     using existing comments in the code or the printed documentation*

//*-*   - The original Minuit subroutines are now member functions.       *

//*-*   - Constructors and destructor have been added.                    *

//*-*   - Instead of passing the FCN  function in the argument            *

//*-*     list, the addresses of this function is stored as pointer       *

//*-*     in the data members of the class. This is by far more elegant   *

//*-*     and flexible in an interactive environment.                     *

//*-*     The member function SetFCN can be used to define this pointer.  *

//*-*   - The ROOT static function Printf is provided to replace all      *

//*-*     format statements and to print on currently defined output file.*

//*-*   - The derived class MinuitOld contains obsolete routines from    *

//*-*     the Fortran based version.                                      *

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



#include <cstdlib>

#include <cstdio>

#include <cstring>

#include <cstdarg>

#include <cmath>



#include "CalRecon/Midnight.h"



#define TMath_Min(a,b) ((a)<(b)?a:b)

#define TMath_Max(a,b) ((a)>(b)?a:b)

#define TMath_Log(a) log(a)

#define TMath_Log10(a) log10(a)

#define TMath_Sqrt(a) sqrt(a)

#define TMath_Power(a,b) pow(a,b)

#define TMath_Abs(a) fabs(a)

#define TMath_Sin(a) sin(a)

#define TMath_Cos(a) cos(a)

#define TMath_ATan(a) atan(a)

#define TMath_ASin(a) asin(a)



static char sBuffer[1024];

static char* Form(char*,...);

static void Printf(const char*,...);



const MInt kMAXDIM = 50;

const MBool kTRUE   = 1;

const MBool kFALSE  = 0;



//______________________________________________________________________________

Midnight::Midnight()

{

//*-*-*-*-*-*-*-*-*-*-*Minuit normal constructor*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

//*-*                  ========================



   fEmpty = 1;

   fPrintf = Printf;



}



//______________________________________________________________________________

Midnight::Midnight(MInt maxpar)

{

//*-*-*-*-*-*-*-*-*-*-*Minuit normal constructor*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

//*-*                  ========================



   BuildArrays(maxpar);



   fEmpty      = 0;

   fPrintf = Printf;



   mninit(5,6,7);

}



//______________________________________________________________________________

Midnight::Midnight(const Midnight &)

{

   // Private copy ctor. This object can not be copied.



  //debug Error("Midnight", "can not copy construct Midnight");

   fPrintf("Midnight can not copy construct Midnight");

}



//______________________________________________________________________________

Midnight::~Midnight()

{

//*-*-*-*-*-*-*-*-*-*-*Minuit default destructor*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

//*-*                  =========================



   DeleteArrays();

}





// added by R Terrier



int Midnight::GetParameter( int parNo, double &currentValue, double &currentError )

{

// return parameter value and error

     int    err;

     MString  name; // ignored

     double bnd1, bnd2; // ignored



     mnpout( parNo, name, currentValue, currentError, bnd1, bnd2, err );



     return err;

}





//______________________________________________________________________________

void Midnight::BuildArrays(MInt maxpar)

{

//*-*-*-*-*-*-*Create internal Minuit arrays for the maxpar parameters*-*-*-*

//*-*          =======================================================



   MInt mni = 25;

   if (maxpar > 10 && maxpar < 200) mni = maxpar;

   fMaxpar      = mni;

   MInt mnihl  = mni*(mni+1)/2;

   MInt maxcpt = 101;

   MInt mne    = 2*mni;

   fCpnam  = new MString[mne];

   fU      = new MDouble[mne];

   fAlim   = new MDouble[mne];

   fBlim   = new MDouble[mne];

   fErp    = new MDouble[mni];

   fErn    = new MDouble[mni];

   fWerr   = new MDouble[mni];

   fGlobcc = new MDouble[mni];

   fNvarl  = new MInt[mne];

   fNiofex = new MInt[mne];

   fNexofi = new MInt[mne];

   fX      = new MDouble[mni];

   fXt     = new MDouble[mni];

   fDirin  = new MDouble[mni];

   fXs     = new MDouble[mni];

   fXts    = new MDouble[mni];

   fDirins = new MDouble[mni];

   fGrd    = new MDouble[mni];

   fG2     = new MDouble[mni];

   fGstep  = new MDouble[mni];

   fGin    = new MDouble[mni];

   fDgrd   = new MDouble[mni];

   fGrds   = new MDouble[mni];

   fG2s    = new MDouble[mni];

   fGsteps = new MDouble[mni];

   fIpfix  = new MInt[mni];

   fVhmat  = new MDouble[mnihl];

   fVthmat = new MDouble[mnihl];

   fP      = new MDouble[mni*(mni+1)];

   fPstar  = new MDouble[mni];

   fPstst  = new MDouble[mni];

   fPbar   = new MDouble[mni];

   fPrho   = new MDouble[mni];

   fWord7  = new MDouble[30];

   fXpt    = new MDouble[maxcpt];

   fYpt    = new MDouble[maxcpt];

   fChpt   = new MString[maxcpt];

   fOrigin = new MString[100];

   fWarmes = new MString[100];



   for (int i = 0; i < fMaxpar; i++) {

      fErp[i] = 0;

      fErn[i] = 0;

   }

}



//______________________________________________________________________________

void Midnight::DeleteArrays()

{

//*-*-*-*-*-*-*-*-*-*-*-*Delete internal Minuit arrays*-*-*-*-*-*-*-*-*

//*-*                    =============================

   if (fEmpty) return;

   delete [] fCpnam;

   delete [] fU;

   delete [] fAlim;

   delete [] fBlim;

   delete [] fErp;

   delete [] fErn;

   delete [] fWerr;

   delete [] fGlobcc;

   delete [] fNvarl;

   delete [] fNiofex;

   delete [] fNexofi;

   delete [] fX;

   delete [] fXt;

   delete [] fDirin;

   delete [] fXs;

   delete [] fXts;

   delete [] fDirins;

   delete [] fGrd;

   delete [] fG2;

   delete [] fGstep;

   delete [] fGin;

   delete [] fDgrd;

   delete [] fGrds;

   delete [] fG2s;

   delete [] fGsteps;

   delete [] fIpfix;

   delete [] fVhmat;

   delete [] fVthmat;

   delete [] fP;

   delete [] fPstar;

   delete [] fPstst;

   delete [] fPbar;

   delete [] fPrho;

   delete [] fWord7;

   delete [] fXpt;

   delete [] fYpt;

   delete [] fChpt;

   delete [] fOrigin;

   delete [] fWarmes;

   fEmpty = 1;

}



//______________________________________________________________________________

void Midnight::SetFCN(MFunction fcn)

{

//*-*-*-*-*-*-*To set the address of the minimization function*-*-*-*-*-*-*-*

//*-*          ===============================================

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

   fFCN = fcn;

}



void Midnight::SetPrintf(MPrintf aPrintf)

{

//*-*-*-*-*-*-*To set the address of the minimization function*-*-*-*-*-*-*-*

//*-*          ===============================================

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

   fPrintf = aPrintf;

}



//______________________________________________________________________________

void Midnight::mnamin()

{

//*-*-*-*-*-*-*-*-*-*-*-*-*Initialize AMIN*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

//*-*                      ===============

//*-*C        Called  from many places.  Initializes the value of AMIN by

//*-*C        calling the user function. Prints out the function value and

//*-*C        parameter values if Print Flag value is high enough.

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    /* Local variables */

    static MDouble fnew;

    static MInt nparx;



    nparx = fNpar;

    if (fISW[4] >= 1) {

	fPrintf(" FIRST CALL TO USER FUNCTION AT NEW START POINT, WITH IFLAG=4.");

    }

    mnexin(fX);

    (*fFCN)(nparx, fGin, fnew, fU, 4);    ++fNfcn;

    fAmin = fnew;

    fEDM  = fBigedm;

} /* mnamin_ */



//______________________________________________________________________________

void Midnight::mnbins(MDouble a1, MDouble a2, MInt naa, MDouble &bl, MDouble &bh, MInt &nb, MDouble &bwid)

{

//*-*-*-*-*-*-*-*-*-*-*Compute reasonable histogram intervals*-*-*-*-*-*-*-*-*

//*-*                  ======================================

//*-*        Function TO DETERMINE REASONABLE HISTOGRAM INTERVALS

//*-*        GIVEN ABSOLUTE UPPER AND LOWER BOUNDS  A1 AND A2

//*-*        AND DESIRED MAXIMUM NUMBER OF BINS NAA

//*-*        PROGRAM MAKES REASONABLE BINNING FROM BL TO BH OF WIDTH BWID

//*-*        F. JAMES,   AUGUST, 1974 , stolen for Minuit, 1988

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    /* Local variables */

    MDouble awid,ah, al, sigfig, sigrnd, alb;

    static MInt kwid, lwid, na, log_;



    al = TMath_Min(a1,a2);

    ah = TMath_Max(a1,a2);

    if (al == ah) ah = al + 1;



//*-*-       IF NAA .EQ. -1 , PROGRAM USES BWID INPUT FROM CALLING ROUTINE

    if (naa == -1) goto L150;

L10:

    na = naa - 1;

    if (na < 1) na = 1;



//*-*-        GET NOMINAL BIN WIDTH IN EXPON FORM

L20:

    awid = (ah-al) / MDouble(na);

    log_ = MInt(TMath_Log10(awid));

    if (awid <= 1) --log_;

    sigfig = awid*TMath_Power(10., -log_);

//*-*-       ROUND MANTISSA UP TO 2, 2.5, 5, OR 10

    if (sigfig > 2) goto L40;

    sigrnd = 2;

    goto L100;

L40:

    if (sigfig > 2.5) goto L50;

    sigrnd = 2.5;

    goto L100;

L50:

    if (sigfig > 5) goto L60;

    sigrnd = 5;

    goto L100;

L60:

    sigrnd = 1;

    ++log_;

L100:

    bwid = sigrnd*TMath_Power(10., log_);

    goto L200;

//*-*-       GET NEW BOUNDS FROM NEW WIDTH BWID

L150:

    if (bwid <= 0) goto L10;

L200:

    alb  = al / bwid;

    lwid = MInt(alb);

    if (alb < 0) --lwid;

    bl   = bwid*MDouble(lwid);

    alb  = ah / bwid + 1;

    kwid = MInt(alb);

    if (alb < 0) --kwid;

    bh = bwid*MDouble(kwid);

    nb = kwid - lwid;

    if (naa > 5) goto L240;

    if (naa == -1) return;

//*-*-        REQUEST FOR ONE BIN IS DIFFICULT CASE

    if (naa > 1 || nb == 1) return;

    bwid *= 2;

    nb = 1;

    return;

L240:

    if (nb << 1 != naa) return;

    ++na;

    goto L20;

} /* mnbins_ */



//______________________________________________________________________________

void Midnight::mncalf(MDouble *pvec, MDouble &ycalf)

{

//*-*-*-*-*-*-*-*-*-*Transform FCN to find further minima*-*-*-*-*-*-*-*-*-*

//*-*                ====================================

//*-*        Called only from MNIMPR.  Transforms the function FCN

//*-*        by dividing out the quadratic part in order to find further

//*-*        minima.    Calculates  ycalf = (f-fmin)/(x-xmin)*v*(x-xmin)

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    /* Local variables */

    static MInt ndex, i, j, m, n, nparx;

    static MDouble denom, f;



    nparx = fNpar;

    mninex(&pvec[0]);

    (*fFCN)(nparx, fGin, f, fU, 4);    ++fNfcn;

    for (i = 1; i <= fNpar; ++i) {

	fGrd[i-1] = 0;

	for (j = 1; j <= fNpar; ++j) {

	    m = TMath_Max(i,j);

	    n = TMath_Min(i,j);

	    ndex = m*(m-1) / 2 + n;

	    fGrd[i-1] += fVthmat[ndex-1]*(fXt[j-1] - pvec[j-1]);

	}

    }

    denom = 0;

    for (i = 1; i <= fNpar; ++i) {denom += fGrd[i-1]*(fXt[i-1] - pvec[i-1]); }

    if (denom <= 0) {

	fDcovar = 1;

	fISW[1] = 0;

	denom   = 1;

    }

    ycalf = (f - fApsi) / denom;

} /* mncalf_ */



//______________________________________________________________________________

void Midnight::mncler()

{

//*-*-*-*-*-*-*-*-*-*-*Resets the parameter list to UNDEFINED*-*-*-*-*-*-*-*

//*-*                  ======================================

//*-*        Called from MINUIT and by option from MNEXCM

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    MInt i;



    fNpfix = 0;

    fNu = 0;

    fNpar = 0;

    fNfcn = 0;

    fNwrmes[0] = 0;

    fNwrmes[1] = 0;

    for (i = 1; i <= fMaxext; ++i) {

	fU[i-1]      = 0;

	fCpnam[i-1]  = fCundef;

	fNvarl[i-1]  = -1;

	fNiofex[i-1] = 0;

    }

    mnrset(1);

    fCfrom  = "CLEAR   ";

    fNfcnfr = fNfcn;

    fCstatu = "UNDEFINED ";

    fLnolim = kTRUE;

    fLphead = kTRUE;

} /* mncler_ */



//______________________________________________________________________________

void Midnight::mncntr(MInt ke1, MInt ke2, MInt &ierrf)

{

//*-*-*-*-*Print function contours in two variables, on line printer*-*-*-*-*

//*-*      =========================================================

//*-*

//*-*                input arguments: parx, pary, devs, ngrid

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    static MString clabel = "0123456789ABCDEFGHIJ";



    /* Local variables */

    MDouble d__1, d__2;

    MDouble fcna[115], fcnb[115], contur[20];

    MDouble  ylabel, fmn, fmx, xlo, ylo, xup, yup;

    MDouble devs, xsav, ysav,  bwidx,  bwidy, unext, ff, xb4;

    static MInt i,  ngrid, ixmid, nparx, ix, nx, ny, ki1, ki2, ixzero, iy, ics;

    static MString chmid, chln, chzero;



    if (ke1 <= 0 || ke2 <= 0) goto L1350;

    if (ke1 > fNu || ke2 > fNu) goto L1350;

    ki1 = fNiofex[ke1-1];

    ki2 = fNiofex[ke2-1];

    if (ki1 <= 0 || ki2 <= 0) goto L1350;

    if (ki1 == ki2) goto L1350;



    if (fISW[1] < 1) {

	mnhess();

	mnwerr();

    }

    nparx = fNpar;

    xsav = fU[ke1-1];

    ysav = fU[ke2-1];

    devs = fWord7[2];

    if (devs <= 0) devs = 2;

    xlo = fU[ke1-1] - devs*fWerr[ki1-1];

    xup = fU[ke1-1] + devs*fWerr[ki1-1];

    ylo = fU[ke2-1] - devs*fWerr[ki2-1];

    yup = fU[ke2-1] + devs*fWerr[ki2-1];

    ngrid = MInt(fWord7[3]);

    if (ngrid <= 0) {

	ngrid = 25;

//*-*  Computing MIN

	nx = TMath_Min(fNpagwd - 15,ngrid);

//*-*  Computing MIN

	ny = TMath_Min(fNpagln - 7,ngrid);

    } else {

	nx = ngrid;

	ny = ngrid;

    }

    if (nx < 11)   nx = 11;

    if (ny < 11)   ny = 11;

    if (nx >= 115) nx = 114;



//*-*-        ask if parameter outside limits

    if (fNvarl[ke1-1] > 1) {

	if (xlo < fAlim[ke1-1]) xlo = fAlim[ke1-1];

	if (xup > fBlim[ke1-1]) xup = fBlim[ke1-1];

    }

    if (fNvarl[ke2-1] > 1) {

	if (ylo < fAlim[ke2-1]) ylo = fAlim[ke2-1];

	if (yup > fBlim[ke2-1]) yup = fBlim[ke2-1];

    }

    bwidx = (xup - xlo) / MDouble(nx);

    bwidy = (yup - ylo) / MDouble(ny);

    ixmid = MInt(((xsav - xlo)*MDouble(nx) / (xup - xlo)) + 1);

    if (fAmin == fUndefi) mnamin();



    for (i = 1; i <= 20; ++i) {	contur[i-1] = fAmin + fUp*(i-1)*(i-1); }

    contur[0] += fUp*.01;

//*-*-               fill FCNB to prepare first row, and find column zero/

    fU[ke2-1] = yup;

    ixzero = 0;

    xb4 = 1;

//TH

    chmid.resize(nx+1);

    chzero.resize(nx+1);

    chln.resize(nx+1);

    for (ix = 1; ix <= nx + 1; ++ix) {

	fU[ke1-1] = xlo + MDouble(ix-1)*bwidx;

	(*fFCN)(nparx, fGin, ff, fU, 4);

	fcnb[ix-1] = ff;

	if (xb4 < 0 && fU[ke1-1] > 0) ixzero = ix - 1;

	xb4          = fU[ke1-1];

	chmid[ix-1]  = '*';

	chzero[ix-1] = '-';

    }

    fPrintf(" Y-AXIS: PARAMETER %3d: %s",ke2,(const char*)fCpnam[ke2-1]);

    if (ixzero > 0) {

	chzero[ixzero-1] = '+';

	chln = " ";

	fPrintf("             X=0");

    }

//*-*-                loop over rows

    for (iy = 1; iy <= ny; ++iy) {

	unext = fU[ke2-1] - bwidy;

//*-*-                prepare this line background pattern for contour

	chln = " ";

// TH	

	chln.resize(nx+1);

	chln[ixmid-1] = '*';

	if (ixzero != 0) chln[ixzero-1] = ':';

	if (fU[ke2-1] > ysav && unext < ysav) chln = chmid;

	if (fU[ke2-1] > 0 && unext < 0)       chln = chzero;

	fU[ke2-1] = unext;

	ylabel = fU[ke2-1] + bwidy*.5;

//*-*-                move FCNB to FCNA and fill FCNB with next row

	for (ix = 1; ix <= nx + 1; ++ix) {

	    fcna[ix-1] = fcnb[ix-1];

	    fU[ke1-1] = xlo + MDouble(ix-1)*bwidx;

	    (*fFCN)(nparx, fGin, ff, fU, 4);

	    fcnb[ix-1] = ff;

	}

//*-*-                look for contours crossing the FCNxy squares

	for (ix = 1; ix <= nx; ++ix) {

	    d__1 = TMath_Max(fcna[ix-1],fcnb[ix-1]),

	    d__2 = TMath_Max(fcna[ix],fcnb[ix]);

	    fmx  = TMath_Max(d__1,d__2);

	    d__1 = TMath_Min(fcna[ix-1],fcnb[ix-1]),

            d__2 = TMath_Min(fcna[ix],fcnb[ix]);

	    fmn  = TMath_Min(d__1,d__2);

	    for (ics = 1; ics <= 20; ++ics) {

		if (contur[ics-1] > fmn)  goto L240;

	    }

	    continue;

L240:

	    if (contur[ics-1] < fmx) chln[ix-1] = clabel[ics-1];

	}

//*-*-                print a row of the contour plot

        fPrintf(" %12.4g %s",ylabel,(const char*)chln);

    }

//*-*-                contours printed, label x-axis

    chln          = " ";

    chln[0]       = 'I';

    chln[ixmid-1] = 'I';

    chln[nx-1]    = 'I';

    fPrintf("              %s",(const char*)chln);



//*-*-               the hardest of all: print x-axis scale!

    chln =  " ";

    if (nx <= 26) {

        fPrintf("        %12.4g%s%12.4g",xlo,(const char*)chln,xup);

        fPrintf("              %s%12.4g",(const char*)chln,xsav);

    } else {

        fPrintf("        %12.4g%s%12.4g%s%12.4g",xlo,(const char*)chln,xsav,(const char*)chln,xup);

    }

    fPrintf("       X-AXIS: PARAMETER%3d: %s  ONE COLUMN=%12.4g"

            ,ke1,(const char*)fCpnam[ke1-1],bwidx);

    fPrintf(" FUNCTION VALUES: F(I)=%12.4g +%12.4g *I**2",fAmin,fUp);

//*-*-                finished.  reset input values

    fU[ke1-1] = xsav;

    fU[ke2-1] = ysav;

    ierrf     = 0;

    return;

L1350:

    fPrintf(" INVALID PARAMETER NUMBER(S) REQUESTED.  IGNORED.");

    ierrf = 1;

} /* mncntr_ */



//______________________________________________________________________________

void Midnight::mncomd(MString crdbin, MInt &icondn)

{

//*-*-*-*-*-*-*-*-*-*-*Reads a command string and executes*-*-*-*-*-*-*-*-*-*

//*-*                  ===================================

//*-*        Called by user.  'Reads' a command string and executes.

//*-*     Equivalent to MNEXCM except that the command is given as a

//*-*          character string.

//*-*

//*-*     ICONDN = 0: command executed normally

//*-*              1: command is blank, ignored

//*-*              2: command line unreadable, ignored

//*-*              3: unknown command, ignored

//*-*              4: abnormal termination (e.g., MIGRAD not converged)

//*-*              5: command is a request to read PARAMETER definitions

//*-*              6: 'SET INPUT' command

//*-*              7: 'SET TITLE' command

//*-*              8: 'SET COVAR' command

//*-*              9: reserved

//*-*             10: END command

//*-*             11: EXIT or STOP command

//*-*             12: RETURN command

//*-*

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    static MString clower = "abcdefghijklmnopqrstuvwxyz";

    static MString cupper = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";



    /* Local variables */

    static MDouble plist[30];

    static MInt ierr, ipos, i, llist, ic, lenbuf, lnc;

    static MBool leader;

    static MString comand, crdbuf, ctemp;



    lenbuf = strlen((const char*)crdbin);

    crdbuf = crdbin;

    icondn = 0;

//*-*-    record not case-sensitive, get upper case, strip leading blanks

    leader = kTRUE;

    ipos = 1;

    for (i = 1; i <= TMath_Min(20,lenbuf); ++i) {

	if (crdbuf[i-1] == '\'') goto L111;

	if (crdbuf[i-1] == ' ') {

	    if (leader) ++ipos;

	    continue;

	}

	leader = kFALSE;

	for (ic = 1; ic <= 26; ++ic) {

	    if (crdbuf[i-1] == clower[ic-1]) crdbuf[i-1] = cupper[ic-1];

	}

    }

L111:

//*-*-                    blank or null command

    if (ipos > lenbuf) {

	fPrintf(" BLANK COMMAND IGNORED.");

	icondn = 1;

	goto L900;

    }

//*-*-                                          . .   preemptive commands

//*-*-              if command is 'PARAMETER'

    if (crdbuf(ipos-1,3) == "PAR") {

	icondn = 5;

	fLphead = kTRUE;

	goto L900;

    }

//*-*-              if command is 'SET INPUT'

    if (crdbuf(ipos-1,3) == "SET INP") {

	icondn = 6;

	fLphead = kTRUE;

	goto L900;

    }

//*-*-              if command is 'SET TITLE'

    if (crdbuf(ipos-1,7) == "SET TIT") {

	icondn = 7;

	fLphead = kTRUE;

	goto L900;

    }

//*-*-              if command is 'SET COVARIANCE'

    if (crdbuf(ipos-1,7) == "SET COV") {

	icondn = 8;

	fLphead = kTRUE;

	goto L900;

    }

//*-*-              crack the command . . . . . . . . . . . . . . . .

    ctemp = crdbuf(ipos-1,7);

    mncrck(ctemp, 20, comand, lnc, 30, plist, llist, ierr, fIsyswr);

    if (ierr > 0) {

	fPrintf(" COMMAND CANNOT BE INTERPRETED");

	icondn = 2;

	goto L900;

    }



    mnexcm(comand, plist, llist, ierr);

    icondn = ierr;

L900:

    return;

} /* mncomd_ */



//______________________________________________________________________________

void Midnight::mncont(MInt ke1, MInt ke2, MInt nptu, MDouble *xptu, MDouble *yptu, MInt &ierrf)

{

//*-*-*-*-*-*-*Find points along a contour where FCN is minimum*-*-*-*-*-*-*

//*-*          ================================================

//*-*       Find NPTU points along a contour where the function

//*-*             FMIN (X(KE1),X(KE2)) =  AMIN+UP

//*-*       where FMIN is the minimum of FCN with respect to all

//*-*       the other NPAR-2 variable parameters (if any).

//*-*   IERRF on return will be equal to the number of points found:

//*-*     NPTU if normal termination with NPTU points found

//*-*     -1   if errors in the calling sequence (KE1, KE2 not variable)

//*-*      0   if less than four points can be found (using MNMNOT)

//*-*     n>3  if only n points can be found (n < NPTU)

//*-*

//*-*                 input arguments: parx, pary, devs, ngrid

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* System generated locals */

    MInt i__1;



    /* Local variables */

    MDouble d__1, d__2;

    static MDouble dist, xdir, ydir, aopt,  u1min, u2min;

    static MDouble gcc[kMAXDIM], w[kMAXDIM], abest, scalx, scaly;

    static MDouble a1, a2, val2mi, val2pl, dc, sclfac, bigdis, sigsav;

    static MInt nall, iold, line, mpar, ierr, inew, move, next, i, j, nfcol, iercr;

    static MInt idist, npcol, kints, i2, i1, lr, nfcnco, ki1, ki2, ki3, ke3;

    static MInt nowpts, istrav, nfmxin, isw2, isw4;

    static MBool ldebug;



    /* Function Body */

    ldebug = fIdbg[6] >= 1;

    if (ke1 <= 0 || ke2 <= 0) goto L1350;

    if (ke1 > fNu || ke2 > fNu) goto L1350;

    ki1 = fNiofex[ke1-1];

    ki2 = fNiofex[ke2-1];

    if (ki1 <= 0 || ki2 <= 0) goto L1350;

    if (ki1 == ki2) goto L1350;

    if (nptu < 4)  goto L1400;



    nfcnco  = fNfcn;

    fNfcnmx = (nptu + 5)*100*(fNpar + 1);

//*-*-          The minimum

    mncuve();

    u1min  = fU[ke1-1];

    u2min  = fU[ke2-1];

    ierrf  = 0;

    fCfrom = "MNContour ";

    fNfcnfr = nfcnco;

    if (fISW[4] >= 0) {

        fPrintf(" START MNCONTOUR CALCULATION OF%4d POINTS ON CONTOUR.",nptu);

	if (fNpar > 2) {

	    if (fNpar == 3) {

		ki3 = 6 - ki1 - ki2;

		ke3 = fNexofi[ki3-1];

		fPrintf(" EACH POINT IS A MINIMUM WITH RESPECT TO PARAMETER %3d  %s",ke3,(const char*)fCpnam[ke3-1]);

	    } else {

		fPrintf(" EACH POINT IS A MINIMUM WITH RESPECT TO THE OTHER%3d VARIABLE PARAMETERS.",fNpar - 2);

	    }

	}

    }



//*-*-          Find the first four points using MNMNOT

//*-*-             ........................ first two points

    mnmnot(ke1, ke2, val2pl, val2mi);

    if (fErn[ki1-1] == fUndefi) {

	xptu[0] = fAlim[ke1-1];

	mnwarn("W", "MNContour ", "Contour squeezed by parameter limits.");

    } else {

	if (fErn[ki1-1] >= 0) goto L1500;

	xptu[0] = u1min + fErn[ki1-1];

    }

    yptu[0] = val2mi;



    if (fErp[ki1-1] == fUndefi) {

	xptu[2] = fBlim[ke1-1];

	mnwarn("W", "MNContour ", "Contour squeezed by parameter limits.");

    } else {

	if (fErp[ki1-1] <= 0) goto L1500;

	xptu[2] = u1min + fErp[ki1-1];

    }

    yptu[2] = val2pl;

    scalx = 1 / (xptu[2] - xptu[0]);

//*-*-             ........................... next two points

    mnmnot(ke2, ke1, val2pl, val2mi);

    if (fErn[ki2-1] == fUndefi) {

	yptu[1] = fAlim[ke2-1];

	mnwarn("W", "MNContour ", "Contour squeezed by parameter limits.");

    } else {

	if (fErn[ki2-1] >= 0) goto L1500;

	yptu[1] = u2min + fErn[ki2-1];

    }

    xptu[1] = val2mi;

    if (fErp[ki2-1] == fUndefi) {

	yptu[3] = fBlim[ke2-1];

	mnwarn("W", "MNContour ", "Contour squeezed by parameter limits.");

    } else {

	if (fErp[ki2-1] <= 0) goto L1500;

	yptu[3] = u2min + fErp[ki2-1];

    }

    xptu[3] = val2pl;

    scaly   = 1 / (yptu[3] - yptu[1]);

    nowpts  = 4;

    next    = 5;

    if (ldebug) {

	fPrintf(" Plot of four points found by MINOS");

	fXpt[0]  = u1min;

	fYpt[0]  = u2min;

	fChpt[0] = ' ';

//*-*  Computing MIN

	nall = TMath_Min(nowpts + 1,101);

	for (i = 2; i <= nall; ++i) {

	    fXpt[i-1] = xptu[i-2];

	    fYpt[i-1] = yptu[i-2];

	}

	fChpt[1] = 'A';

	fChpt[2] = 'B';

	fChpt[3] = 'C';

	fChpt[4] = 'D';

	mnplot(fXpt, fYpt, fChpt, nall, fNpagwd, fNpagln);

    }



//*-*-              ..................... save some values before fixing

    isw2   = fISW[1];

    isw4   = fISW[3];

    sigsav = fEDM;

    istrav = fIstrat;

    dc     = fDcovar;

    fApsi  = fEpsi*.5;

    abest  = fAmin;

    mpar   = fNpar;

    nfmxin = fNfcnmx;

    for (i = 1; i <= mpar; ++i) { fXt[i-1] = fX[i-1]; }

    i__1 = mpar*(mpar + 1) / 2;

    for (j = 1; j <= i__1; ++j) { fVthmat[j-1] = fVhmat[j-1]; }

    for (i = 1; i <= mpar; ++i) {

	gcc[i-1] = fGlobcc[i-1];

	w[i-1] = fWerr[i-1];

    }

//*-*-                          fix the two parameters in question

    kints = fNiofex[ke1-1];

    mnfixp(kints-1, ierr);

    kints = fNiofex[ke2-1];

    mnfixp(kints-1, ierr);

//*-*-              ......................Fill in the rest of the points

    for (inew = next; inew <= nptu; ++inew) {

//*-*            find the two neighbouring points with largest separation

	bigdis = 0;

	for (iold = 1; iold <= inew - 1; ++iold) {

	    i2 = iold + 1;

	    if (i2 == inew) i2 = 1;

	    d__1 = scalx*(xptu[iold-1] - xptu[i2-1]);

	    d__2 = scaly*(yptu[iold-1] - yptu[i2-1]);

	    dist = d__1*d__1 + d__2*d__2;

	    if (dist > bigdis) {

		bigdis = dist;

		idist  = iold;

	    }

	}

	i1 = idist;

	i2 = i1 + 1;

	if (i2 == inew) i2 = 1;

//*-*-                  next point goes between I1 and I2

	a1 = .5;

	a2 = .5;

L300:

	fXmidcr = a1*xptu[i1-1] + a2*xptu[i2-1];

	fYmidcr = a1*yptu[i1-1] + a2*yptu[i2-1];

	xdir    = yptu[i2-1] - yptu[i1-1];

	ydir    = xptu[i1-1] - xptu[i2-1];

	sclfac  = TMath_Max(TMath_Abs(xdir*scalx),TMath_Abs(ydir*scaly));

	fXdircr = xdir / sclfac;

	fYdircr = ydir / sclfac;

	fKe1cr  = ke1;

	fKe2cr  = ke2;

//*-*-               Find the contour crossing point along DIR

	fAmin = abest;

	mncros(aopt, iercr);

	if (iercr > 1) {

//*-*-             If cannot find mid-point, try closer to point 1

	    if (a1 > .5) {

		if (fISW[4] >= 0) {

		    fPrintf(" MNCONT CANNOT FIND NEXT POINT ON CONTOUR.  ONLY%3d POINTS FOUND.",nowpts);

		}

		goto L950;

	    }

	    mnwarn("W", "MNContour ", "Cannot find midpoint, try closer.");

	    a1 = .75;

	    a2 = .25;

	    goto L300;

	}

//*-*-               Contour has been located, insert new point in list

	for (move = nowpts; move >= i1 + 1; --move) {

	    xptu[move] = xptu[move-1];

	    yptu[move] = yptu[move-1];

	}

	++nowpts;

	xptu[i1] = fXmidcr + fXdircr*aopt;

	yptu[i1] = fYmidcr + fYdircr*aopt;

    }

L950:



    ierrf = nowpts;

    fCstatu = "SUCCESSFUL";

    if (nowpts < nptu) 	fCstatu = "INCOMPLETE";



//*-*-               make a lineprinter plot of the contour

    if (fISW[4] >= 0) {

	fXpt[0]  = u1min;

	fYpt[0]  = u2min;

	fChpt[0] = ' ';

	nall = TMath_Min(nowpts + 1,101);

	for (i = 2; i <= nall; ++i) {

	    fXpt[i-1]  = xptu[i-2];

	    fYpt[i-1]  = yptu[i-2];

	    fChpt[i-1] = 'X';

	}

        fPrintf(" Y-AXIS: PARAMETER %3d  %s",ke2,(const char*)fCpnam[ke2-1]);



	mnplot(fXpt, fYpt, fChpt, nall, fNpagwd, fNpagln);



	fPrintf("                         X-AXIS: PARAMETER %3d  %s",ke1,(const char*)fCpnam[ke1-1]);

    }

//*-*-                print out the coordinates around the contour

    if (fISW[4] >= 1) {

	npcol = (nowpts + 1) / 2;

	nfcol = nowpts / 2;

	fPrintf("%5d POINTS ON CONTOUR.   FMIN=%13.5e   ERRDEF=%11.3g",nowpts,abest,fUp);

        fPrintf("         %s%s%s%s",(const char*)fCpnam[ke1-1],

                                  (const char*)fCpnam[ke2-1],

                                  (const char*)fCpnam[ke1-1],

                                  (const char*)fCpnam[ke2-1]);

	for (line = 1; line <= nfcol; ++line) {

	    lr = line + npcol;

	    fPrintf(" %5d%13.5e%13.5e          %5d%13.5e%13.5e",line,xptu[line-1],yptu[line-1],lr,xptu[lr-1],yptu[lr-1]);

	}

	if (nfcol < npcol) {

	    fPrintf(" %5d%13.5e%13.5e",npcol,xptu[npcol-1],yptu[npcol-1]);

	}

    }

//*-*-                                   . . contour finished. reset v

    fItaur = 1;

    mnfree(1);

    mnfree(1);

    i__1 = mpar*(mpar + 1) / 2;

    for (j = 1; j <= i__1; ++j) { fVhmat[j-1] = fVthmat[j-1]; }

    for (i = 1; i <= mpar; ++i) {

	fGlobcc[i-1] = gcc[i-1];

	fWerr[i-1]   = w[i-1];

	fX[i-1]      = fXt[i-1];

    }

    mninex(fX);

    fEDM    = sigsav;

    fAmin   = abest;

    fISW[1] = isw2;

    fISW[3] = isw4;

    fDcovar = dc;

    fItaur  = 0;

    fNfcnmx = nfmxin;

    fIstrat = istrav;

    fU[ke1-1] = u1min;

    fU[ke2-1] = u2min;

    goto L2000;

//*-*-                                    Error returns

L1350:

    fPrintf(" INVALID PARAMETER NUMBERS.");

    goto L1450;

L1400:

    fPrintf(" LESS THAN FOUR POINTS REQUESTED.");

L1450:

    ierrf   = -1;

    fCstatu = "USER ERROR";

    goto L2000;

L1500:

    fPrintf(" MNCONT UNABLE TO FIND FOUR POINTS.");

    fU[ke1-1] = u1min;

    fU[ke2-1] = u2min;

    ierrf     = 0;

    fCstatu   = "FAILED";

L2000:

    fCfrom  = "MNContour ";

    fNfcnfr = nfcnco;

} /* mncont_ */



//______________________________________________________________________________

void Midnight::mncrck(MString &crdbuf, MInt maxcwd, MString &comand, MInt &lnc,

        MInt mxp, MDouble *plist, MInt &llist, MInt &ierr, MInt)

{

//*-*-*-*-*-*-*-*-*-*-*-*Cracks the free-format input*-*-*-*-*-*-*-*-*-*-*-*-*

//*-*                    ============================

//*-*       Called from MNREAD.

//*-*       Cracks the free-format input, expecting zero or more

//*-*         alphanumeric fields (which it joins into COMAND(1:LNC))

//*-*         followed by one or more numeric fields separated by

//*-*         blanks and/or one comma.  The numeric fields are put into

//*-*         the LLIST (but at most MXP) elements of PLIST.

//*-*      IERR = 0 if no errors,

//*-*           = 1 if error(s).

//*-*      Diagnostic messages are written to ISYSWR

//*-*

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Initialized data */



    static MString cnull  = ")NULL STRING   ";

    static MString cnumer = "123456789-.0+";



    /* Local variables */

    static MInt ifld, iend, lend, left, nreq, ipos, kcmnd, nextb, ic, ibegin, ltoadd;

    static MInt ielmnt, lelmnt[25], nelmnt;

    static MString celmnt[25], ctemp;



    /* Function Body */

    ielmnt = 0;

    lend   = strlen((const char*)crdbuf);

    nextb  = 1;

    ierr   = 0;

//*-*-                                  . . . .  loop over words CELMNT

L10:

    for (ipos = nextb; ipos <= lend; ++ipos) {

	ibegin = ipos;

	if (crdbuf[ipos-1] == ' ') continue;

	if (crdbuf[ipos-1] == ',') goto L250;

	goto L150;

    }

    goto L300;

L150:

//*-*-              found beginning of word, look for end

    for (ipos = ibegin + 1; ipos <= lend; ++ipos) {

	if (crdbuf[ipos-1] == ' ') goto L250;

	if (crdbuf[ipos-1] == ',') goto L250;

    }

    ipos = lend + 1;

L250:

    iend = ipos - 1;

    ++ielmnt;

    if (iend >= ibegin) celmnt[ielmnt-1] = crdbuf(ibegin-1, iend-ibegin+1);

    else                celmnt[ielmnt-1] = cnull;

    lelmnt[ielmnt-1] = iend - ibegin + 1;

    if (lelmnt[ielmnt-1] > 19) {

	fPrintf(" MINUIT WARNING: INPUT DATA WORD TOO LONG.");

        ctemp = crdbuf(ibegin-1,iend-ibegin+1);

        fPrintf("     ORIGINAL:%s",(const char*)ctemp);

        fPrintf(" TRUNCATED TO:%s",(const char*)celmnt[ielmnt-1]);

	lelmnt[ielmnt-1] = 19;

    }

    if (ipos >= lend) goto L300;

    if (ielmnt >= 25) goto L300;

//*-*-                    look for comma or beginning of next word

    for (ipos = iend + 1; ipos <= lend; ++ipos) {

	if (crdbuf[ipos-1] == ' ') continue;

	nextb = ipos;

	if (crdbuf[ipos-1] == ',') nextb = ipos + 1;

	goto L10;

    }

//*-*-                All elements found, join the alphabetic ones to

//*-*-                               form a command

L300:

    nelmnt   = ielmnt;

    comand   = " ";

    lnc      = 1;

    plist[0] = 0;

    llist    = 0;

    if (ielmnt == 0) goto L900;

    kcmnd = 0;

    for (ielmnt = 1; ielmnt <= nelmnt; ++ielmnt) {

	if ( celmnt[ielmnt-1] == cnull) goto L450;

	for (ic = 1; ic <= 13; ++ic) {

	    if (celmnt[ielmnt-1](0,1) == cnumer(ic-1,1)) goto L450;

	}

	if (kcmnd >= maxcwd) continue;

	left = maxcwd - kcmnd;

	ltoadd = lelmnt[ielmnt-1];

	if (ltoadd > left) ltoadd = left;

	//comand(kcmnd,ltoadd) =  celmnt[ielmnt-1](0,ltoadd);

	comand.replace(kcmnd,ltoadd,celmnt[ielmnt-1](0,ltoadd));

	kcmnd += ltoadd;

	if (kcmnd == maxcwd) continue;

	comand[kcmnd] = ' ';

	++kcmnd;

    }

    lnc = kcmnd;

    goto L900;

L450:

    lnc = kcmnd;

//*-*-                     . . . .  we have come to a numeric field

    llist = 0;

    for (ifld = ielmnt; ifld <= nelmnt; ++ifld) {

	++(llist);

	if (llist > mxp) {

	    nreq = nelmnt - ielmnt + 1;

	    fPrintf(" MINUIT WARNING IN MNCRCK: ");

            fPrintf(" COMMAND HAS INPUT%5d NUMERIC FIELDS, BUT MINUIT CAN ACCEPT ONLY%3d",nreq,mxp);

	    goto L900;

	}

	if (celmnt[ifld-1] == cnull) plist[llist-1] = 0;

	else {

            fPrintf("Fatal Error: mnerr code not yet implemented.");

	}

    }

//*-*-                                 end loop over numeric fields

L900:

    if (lnc <= 0) lnc = 1;

} /* mncrck_ */



//______________________________________________________________________________

void Midnight::mncros(MDouble &aopt, MInt &iercr)

{

//*-*-*-*-*-*-*-*-*-*-*Find point where MNEVAL=AMIN+UP*-*-*-*-*-*-*-*-*-*-*-*

//*-*                  ===============================

//*-*       Find point where MNEVAL=AMIN+UP, along the line through

//*-*       XMIDCR,YMIDCR with direction XDIRCR,YDIRCR,   where X and Y

//*-*       are parameters KE1CR and KE2CR.  If KE2CR=0 (from MINOS),

//*-*       only KE1CR is varied.  From MNCONT, both are varied.

//*-*       Crossing point is at

//*-*        (U(KE1),U(KE2)) = (XMID,YMID) + AOPT*(XDIR,YDIR)

//*-*

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    static MString charal = " .ABCDEFGHIJKLMNOPQRSTUVWXYZ";



    /* Local variables */

    static MDouble alsb[3], flsb[3], bmin, bmax, zmid, sdev, zdir, zlim;

    static MDouble coeff[3], aleft, aulim, fdist, adist, aminsv;

    static MDouble anext, fnext, slope, s1, s2, x1, x2, ecarmn, ecarmx;

    static MDouble determ, rt, smalla, aright, aim, tla, tlf, dfda,ecart;

    static MInt iout, i, ileft, ierev, maxlk, ibest, ik, it;

    static MInt noless, iworst, iright, itoohi, kex, ipt;

    static MBool ldebug;

    static const char *chsign;

    x2 = 0;



    ldebug = fIdbg[6] >= 1;

    aminsv = fAmin;

//*-*-       convergence when F is within TLF of AIM and next prediction

//*-*-       of AOPT is within TLA of previous value of AOPT

    aim      = fAmin + fUp;

    tlf      = fUp*.01;

    tla      = .01;

    fXpt[0]  = 0;

    fYpt[0]  = aim;

    fChpt[0] = ' ';

    ipt = 1;

    if (fKe2cr == 0) {

	fXpt[1]  = -1;

	fYpt[1]  = fAmin;

	fChpt[1] = '.';

	ipt      = 2;

    }

//*-*-                   find the largest allowed A

    aulim = 100;

    for (ik = 1; ik <= 2; ++ik) {

	if (ik == 1) {

	    kex  = fKe1cr;

	    zmid = fXmidcr;

	    zdir = fXdircr;

	} else {

	    if (fKe2cr == 0) continue;

	    kex  = fKe2cr;

	    zmid = fYmidcr;

	    zdir = fYdircr;

	}

	if (fNvarl[kex-1] <= 1) continue;

	if (zdir == 0) continue;

	zlim = fAlim[kex-1];

	if (zdir > 0) zlim = fBlim[kex-1];

	aulim = TMath_Min(aulim,(zlim - zmid) / zdir);

    }

//*-*-                 LSB = Line Search Buffer

//*-*-         first point

    anext   = 0;

    aopt    = anext;

    fLimset = kFALSE;

    if (aulim < aopt + tla) fLimset = kTRUE;

    mneval(anext, fnext, ierev);

//*-* debug printout:

    if (ldebug) {

	fPrintf(" MNCROS: calls=%8d   AIM=%10.5f  F,A=%10.5f%10.5f",fNfcn,aim,fnext,aopt);

    }

    if (ierev > 0) goto L900;

    if (fLimset && fnext <= aim) goto L930;

    ++ipt;

    fXpt[ipt-1]  = anext;

    fYpt[ipt-1]  = fnext;

    fChpt[ipt-1] = charal(ipt-1,1);

    alsb[0] = anext;

    flsb[0] = fnext;

    fnext   = TMath_Max(fnext,aminsv + fUp*.1);

    aopt    = TMath_Sqrt(fUp / (fnext - aminsv)) - 1;

    if (TMath_Abs(fnext - aim) < tlf) goto L800;



    if (aopt < -.5)aopt = -.5;

    if (aopt > 1)  aopt = 1;

    fLimset = kFALSE;

    if (aopt > aulim) {

	aopt    = aulim;

	fLimset = kTRUE;

    }

    mneval(aopt, fnext, ierev);

//*-* debug printout:

    if (ldebug) {

	fPrintf(" MNCROS: calls=%8d   AIM=%10.5f  F,A=%10.5f%10.5f",fNfcn,aim,fnext,aopt);

    }

    if (ierev > 0) goto L900;

    if (fLimset && fnext <= aim) goto L930;

    alsb[1] = aopt;

    ++ipt;

    fXpt[ipt-1]  = alsb[1];

    fYpt[ipt-1]  = fnext;

    fChpt[ipt-1] = charal(ipt-1,1);

    flsb[1]      = fnext;

    dfda         = (flsb[1] - flsb[0]) / (alsb[1] - alsb[0]);

//*-*-                  DFDA must be positive on the contour

    if (dfda > 0) goto L460;

L300:

    mnwarn("D", "MNCROS    ", "Looking for slope of the right sign");

    maxlk = 15 - ipt;

    for (it = 1; it <= maxlk; ++it) {

	alsb[0] = alsb[1];

	flsb[0] = flsb[1];

	aopt    = alsb[0] + MDouble(it)*.2;

	fLimset = kFALSE;

	if (aopt > aulim) {

	    aopt    = aulim;

	    fLimset = kTRUE;

	}

	mneval(aopt, fnext, ierev);

//*-* debug printout:

	if (ldebug) {

            fPrintf(" MNCROS: calls=%8d   AIM=%10.5f  F,A=%10.5f%10.5f",fNfcn,aim,fnext,aopt);

	}

	if (ierev > 0) goto L900;

	if (fLimset && fnext <= aim) goto L930;

	alsb[1] = aopt;

	++ipt;

	fXpt[ipt-1]  = alsb[1];

	fYpt[ipt-1]  = fnext;

	fChpt[ipt-1] = charal(ipt-1,1);

	flsb[1]      = fnext;

	dfda         = (flsb[1] - flsb[0]) / (alsb[1] - alsb[0]);

	if (dfda > 0) goto L450;

    }

    mnwarn("W", "MNCROS    ", "Cannot find slope of the right sign");

    goto L950;

L450:

//*-*-                   we have two points with the right slope

L460:

    aopt  = alsb[1] + (aim - flsb[1]) / dfda;

    fdist = TMath_Min(TMath_Abs(aim - flsb[0]),TMath_Abs(aim - flsb[1]));

    adist = TMath_Min(TMath_Abs(aopt - alsb[0]),TMath_Abs(aopt - alsb[1]));

    tla = .01;

    if (TMath_Abs(aopt) > 1) tla = TMath_Abs(aopt)*.01;

    if (adist < tla && fdist < tlf) goto L800;

    if (ipt >= 15) goto L950;

    bmin = TMath_Min(alsb[0],alsb[1]) - 1;

    if (aopt < bmin) aopt = bmin;

    bmax = TMath_Max(alsb[0],alsb[1]) + 1;

    if (aopt > bmax) aopt = bmax;

//*-*-                   Try a third point

    fLimset = kFALSE;

    if (aopt > aulim) {

	aopt    = aulim;

	fLimset = kTRUE;

    }

    mneval(aopt, fnext, ierev);

//*-* debug printout:

    if (ldebug) {

        fPrintf(" MNCROS: calls=%8d   AIM=%10.5f  F,A=%10.5f%10.5f",fNfcn,aim,fnext,aopt);

    }

    if (ierev > 0) goto L900;

    if (fLimset && fnext <= aim) goto L930;

    alsb[2] = aopt;

    ++ipt;

    fXpt[ipt-1]  = alsb[2];

    fYpt[ipt-1]  = fnext;

    fChpt[ipt-1] = charal(ipt-1,1);

    flsb[2]      = fnext;

//*-*-               now we have three points, ask how many <AIM

    ecarmn = TMath_Abs(fnext-aim);

    ibest  = 3;

    ecarmx = 0;

    noless = 0;

    for (i = 1; i <= 3; ++i) {

	ecart = TMath_Abs(flsb[i-1] - aim);

	if (ecart > ecarmx) { ecarmx = ecart; iworst = i; }

	if (ecart < ecarmn) { ecarmn = ecart; ibest = i; }

	if (flsb[i-1] < aim) ++noless;

    }

//*-*-          if at least one on each side of AIM, fit a parabola

    if (noless == 1 || noless == 2) goto L500;

//*-*-          if all three are above AIM, third must be closest to AIM

    if (noless == 0 && ibest != 3) goto L950;

//*-*-          if all three below, and third is not best, then slope

//*-*-            has again gone negative, look for positive slope.

    if (noless == 3 && ibest != 3) {

	alsb[1] = alsb[2];

	flsb[1] = flsb[2];

	goto L300;

    }

//*-*-          in other cases, new straight line thru last two points

    alsb[iworst-1] = alsb[2];

    flsb[iworst-1] = flsb[2];

    dfda = (flsb[1] - flsb[0]) / (alsb[1] - alsb[0]);

    goto L460;

//*-*-               parabola fit

L500:

    mnpfit(alsb, flsb, 3, coeff, sdev);

    if (coeff[2] <= 0) {

	mnwarn("D", "MNCROS    ", "Curvature is negative near contour line.");

    }

    determ = coeff[1]*coeff[1] - coeff[2]*4*(coeff[0] - aim);

    if (determ <= 0) {

	mnwarn("D", "MNCROS    ", "Problem 2, impossible determinant");

	goto L950;

    }

//*-*-               Find which root is the right one

    rt = TMath_Sqrt(determ);

    x1 = (-coeff[1] + rt) / (coeff[2]*2);

    x2 = (-coeff[1] - rt) / (coeff[2]*2);

    s1 = coeff[1] + x1*2*coeff[2];

    s2 = coeff[1] + x2*2*coeff[2];

    if (s1*s2 > 0) {

	fPrintf(" MNCONTour problem 1");

    }

    aopt  = x1;

    slope = s1;

    if (s2 > 0) {

	aopt  = x2;

	slope = s2;

    }

//*-*-        ask if converged

    tla = .01;

    if (TMath_Abs(aopt) > 1) tla = TMath_Abs(aopt)*.01;

    if (TMath_Abs(aopt - alsb[ibest-1]) < tla && TMath_Abs(flsb[ibest-1] - aim) < tlf) {

	goto L800;

    }

    if (ipt >= 15) goto L950;



//*-*-        see if proposed point is in acceptable zone between L and R

//*-*-        first find ILEFT, IRIGHT, IOUT and IBEST

    ileft  = 0;

    iright = 0;

    ibest  = 1;

    ecarmx = 0;

    ecarmn = TMath_Abs(aim - flsb[0]);

    for (i = 1; i <= 3; ++i) {

	ecart = TMath_Abs(flsb[i-1] - aim);

	if (ecart < ecarmn) { ecarmn = ecart; ibest = i; }

	if (ecart > ecarmx) { ecarmx = ecart; }

	if (flsb[i-1] > aim) {

	    if (iright == 0) iright = i;

	    else if (flsb[i-1] > flsb[iright-1]) iout = i;

	    else { iout = iright; iright = i; }

	}

        else if (ileft == 0) ileft = i;

	else if (flsb[i-1] < flsb[ileft-1]) iout = i;

	else { iout = ileft; ileft = i;	}

    }

//*-*-      avoid keeping a very bad point next time around

    if (ecarmx > TMath_Abs(flsb[iout-1] - aim)*10) {

	aopt = aopt*.5 + (alsb[iright-1] + alsb[ileft-1])*.25;

    }

//*-*-        knowing ILEFT and IRIGHT, get acceptable window

    smalla = tla*.1;

    if (slope*smalla > tlf) smalla = tlf / slope;

    aleft  = alsb[ileft-1] + smalla;

    aright = alsb[iright-1] - smalla;

//*-*-        move proposed point AOPT into window if necessary

    if (aopt < aleft)   aopt = aleft;

    if (aopt > aright)  aopt = aright;

    if (aleft > aright) aopt = (aleft + aright)*.5;



//*-*-        see if proposed point outside limits (should be impossible!)

    fLimset = kFALSE;

    if (aopt > aulim) {

	aopt    = aulim;

	fLimset = kTRUE;

    }

//*-*-                 Evaluate function at new point AOPT

    mneval(aopt, fnext, ierev);

//*-* debug printout:

    if (ldebug) {

        fPrintf(" MNCROS: calls=%8d   AIM=%10.5f  F,A=%10.5f%10.5f",fNfcn,aim,fnext,aopt);

    }

    if (ierev > 0) goto L900;

    if (fLimset && fnext <= aim) goto L930;

    ++ipt;

    fXpt[ipt-1]  = aopt;

    fYpt[ipt-1]  = fnext;

    fChpt[ipt-1] = charal(ipt-1,1);

//*-*-               Replace odd point by new one

    alsb[iout-1] = aopt;

    flsb[iout-1] = fnext;

//*-*-         the new point may not be the best, but it is the only one

//*-*-         which could be good enough to pass convergence criteria

    ibest = iout;

    goto L500;



//*-*-      Contour has been located, return point to MNCONT OR MINOS

L800:

    iercr = 0;

    goto L1000;

//*-*-               error in the minimization

L900:

    if (ierev == 1) goto L940;

    goto L950;

//*-*-               parameter up against limit

L930:

    iercr = 1;

    goto L1000;

//*-*-               too many calls to FCN

L940:

    iercr = 2;

    goto L1000;

//*-*-               cannot find next point

L950:

    iercr = 3;

//*-*-               in any case

L1000:

    if (ldebug) {

	itoohi = 0;

	for (i = 1; i <= ipt; ++i) {

	    if (fYpt[i-1] > aim + fUp) {

		fYpt[i-1]  = aim + fUp;

		fChpt[i-1] = '+';

		itoohi     = 1;

	    }

	}

	chsign = "POSI";

	if (fXdircr < 0) chsign = "NEGA";

	if (fKe2cr == 0) {

	    fPrintf("  %sTIVE MINOS ERROR, PARAMETER %3d",chsign,fKe1cr);

	}

	if (itoohi == 1) {

	    fPrintf("POINTS LABELLED '+' WERE TOO HIGH TO PLOT.");

	}

	if (iercr == 1) {

	    fPrintf("RIGHTMOST POINT IS UP AGAINST LIMIT.");

	}

	mnplot(fXpt, fYpt, fChpt, ipt, fNpagwd, fNpagln);

    }

} /* mncros_ */



//______________________________________________________________________________

void Midnight::mncuve()

{

//*-*-*-*-*-*-*-*Makes sure that the current point is a local minimum*-*-*-*-*

//*-*            ====================================================

//*-*        Makes sure that the current point is a local

//*-*        minimum and that the error matrix exists,

//*-*        or at least something good enough for MINOS and MNCONT

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    /* Local variables */

    static MDouble dxdi, wint;

    static MInt ndex, iext, i, j;



    if (fISW[3] < 1) {

	fPrintf(" FUNCTION MUST BE MINIMIZED BEFORE CALLING %s",(const char*)fCfrom);

	fApsi = fEpsi;

	mnmigr();

    }

    if (fISW[1] < 3) {

	mnhess();

	if (fISW[1] < 1) {

	    mnwarn("W", fCfrom, "NO ERROR MATRIX.  WILL IMPROVISE.");

	    for (i = 1; i <= fNpar; ++i) {

		ndex = i*(i-1) / 2;

		for (j = 1; j <= i-1; ++j) {

		    ++ndex;

		    fVhmat[ndex-1] = 0;

		}

		++ndex;

		if (fG2[i-1] <= 0) {

		    wint = fWerr[i-1];

		    iext = fNexofi[i-1];

		    if (fNvarl[iext-1] > 1) {

			mndxdi(fX[i-1], i-1, dxdi);

			if (TMath_Abs(dxdi) < .001) wint = .01;

			else                   wint /= TMath_Abs(dxdi);

		    }

		    fG2[i-1] = fUp / (wint*wint);

		}

		fVhmat[ndex-1] = 2 / fG2[i-1];

	    }

	    fISW[1] = 1;

	    fDcovar = 1;

	} else  mnwerr();

    }

} /* mncuve_ */



//______________________________________________________________________________

void Midnight::mnderi()

{

//*-*-*-*-*-*-*-*Calculates the first derivatives of FCN (GRD)*-*-*-*-*-*-*-*

//*-*            =============================================

//*-*        Calculates the first derivatives of FCN (GRD),

//*-*        either by finite differences or by transforming the user-

//*-*        supplied derivatives to internal coordinates,

//*-*        according to whether ISW(3) is zero or one.

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    /* Local variables */

    static MDouble step, dfmin, stepb4, dd, df, fs1;

    static MDouble tlrstp, tlrgrd, epspri, optstp, stpmax, stpmin, fs2, grbfor, d1d2, xtf;

    static MInt icyc, ncyc, iint, iext, i, nparx;

    static MBool ldebug;



    nparx = fNpar;

    ldebug = fIdbg[2] >= 1;

    if (fAmin == fUndefi) mnamin();

    if (fISW[2] == 1) goto L100;



    if (ldebug) {

//*-*-                      make sure starting at the right place

	mninex(fX);

	nparx = fNpar;

	(*fFCN)(nparx, fGin, fs1, fU, 4);	++fNfcn;

	if (fs1 != fAmin) {

	    df    = fAmin - fs1;

	    mnwarn("D", "MNDERI", Form("function value differs from AMIN by %12.3g",df));

	    fAmin = fs1;

	}

	fPrintf("  FIRST DERIVATIVE DEBUG PRINTOUT.  MNDERI");

        fPrintf(" PAR    DERIV     STEP      MINSTEP   OPTSTEP  D1-D2    2ND DRV");

    }

    dfmin = fEpsma2*8*(TMath_Abs(fAmin) + fUp);

    if (fIstrat <= 0) {

	ncyc   = 2;

	tlrstp = .5;

	tlrgrd = .1;

    } else if (fIstrat == 1) {

	ncyc   = 3;

	tlrstp = .3;

	tlrgrd = .05;

    } else {

	ncyc   = 5;

	tlrstp = .1;

	tlrgrd = .02;

    }

//*-*-                               loop over variable parameters

    for (i = 1; i <= fNpar; ++i) {

	epspri = fEpsma2 + TMath_Abs(fGrd[i-1]*fEpsma2);

//*-*-        two-point derivatives always assumed necessary

//*-*-        maximum number of cycles over step size depends on strategy

	xtf = fX[i-1];

	stepb4 = 0;

//*-*-                              loop as little as possible here!/

	for (icyc = 1; icyc <= ncyc; ++icyc) {

//*-*-                ........ theoretically best step

	    optstp = TMath_Sqrt(dfmin / (TMath_Abs(fG2[i-1]) + epspri));

//*-*-                    step cannot decrease by more than a factor of ten

	    step = TMath_Max(optstp,TMath_Abs(fGstep[i-1]*.1));

//*-*-                but if parameter has limits, max step size = 0.5

	    if (fGstep[i-1] < 0 && step > .5) step = .5;

//*-*-                and not more than ten times the previous step

	    stpmax = TMath_Abs(fGstep[i-1])*10;

	    if (step > stpmax) step = stpmax;

//*-*-                minimum step size allowed by machine precision

	    stpmin = TMath_Abs(fEpsma2*fX[i-1])*8;

	    if (step < stpmin) step = stpmin;

//*-*-                end of iterations if step change less than factor 2

	    if (TMath_Abs((step - stepb4) / step) < tlrstp) goto L50;

//*-*-        take step positive

	    stepb4 = step;

	    if (fGstep[i-1] > 0) fGstep[i-1] =  TMath_Abs(step);

	    else                 fGstep[i-1] = -TMath_Abs(step);

	    stepb4  = step;

	    fX[i-1] = xtf + step;

	    mninex(fX);

	    (*fFCN)(nparx, fGin, fs1, fU, 4);	    ++fNfcn;

//*-*-        take step negative

	    fX[i-1] = xtf - step;

	    mninex(fX);

	    (*fFCN)(nparx, fGin, fs2, fU, 4);	    ++fNfcn;

	    grbfor = fGrd[i-1];

	    fGrd[i-1] = (fs1 - fs2) / (step*2);

	    fG2[i-1]  = (fs1 + fs2 - fAmin*2) / (step*step);

	    fX[i-1]   = xtf;

	    if (ldebug) {

		d1d2 = (fs1 + fs2 - fAmin*2) / step;

		fPrintf("%4d%11.3g%11.3g%10.2g%10.2g%10.2g%10.2g%10.2g",i,fGrd[i-1],step,stpmin,optstp,d1d2,fG2[i-1]);

	    }

//*-*-        see if another iteration is necessary

	    if (TMath_Abs(grbfor - fGrd[i-1]) / (TMath_Abs(fGrd[i-1]) + dfmin/step) < tlrgrd)

		goto L50;

	}

//*-*-                          end of ICYC loop. too many iterations

	if (ncyc == 1) goto L50;

	mnwarn("D", "MNDERI", Form("First derivative not converged. %g%g",fGrd[i-1],grbfor));

L50:

	;

    }

    mninex(fX);

    return;

//*-*-                                       .  derivatives calc by fcn

L100:

    for (iint = 1; iint <= fNpar; ++iint) {

	iext = fNexofi[iint-1];

	if (fNvarl[iext-1] <= 1) {

           fGrd[iint-1] = fGin[iext-1];

	} else {

           dd = (fBlim[iext-1] - fAlim[iext-1])*.5*TMath_Cos(fX[iint-1]);

           fGrd[iint-1] = fGin[iext-1]*dd;

        }

    }

} /* mnderi_ */



//______________________________________________________________________________

void Midnight::mndxdi(MDouble pint, MInt ipar, MDouble &dxdi)

{

//*-*-*-*Calculates the transformation factor between ext/internal values*-*

//*-*    =====================================================================

//*-*        calculates the transformation factor between external and

//*-*        internal parameter values.     this factor is one for

//*-*        parameters which are not limited.     called from MNEMAT.

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    MInt i = fNexofi[ipar];

    dxdi = 1;

    if (fNvarl[i-1] > 1) {

	dxdi = TMath_Abs((fBlim[i-1] - fAlim[i-1])*TMath_Cos(pint))*.5;

    }

} /* mndxdi_ */



//______________________________________________________________________________

void Midnight::mneig(MDouble *a, MInt ndima, MInt n, MInt mits, MDouble *work, MDouble precis, MInt &ifault)

{

//*-*-*-*-*-*-*-*-*-*-*-*Compute matrix eigen values*-*-*-*-*-*-*-*-*-*-*-*-*

//*-*                    ===========================

    /* System generated locals */

    MInt a_offset;

    MDouble d__1;



    /* Local variables */

    static MDouble b, c, f, h, r, s, hh, gl, pr, pt;

    static MInt i, j, k, l, m, i0, i1, j1, m1, n1;



//*-*-         PRECIS is the machine precision EPSMAC

    /* Parameter adjustments */

    a_offset = ndima + 1;

    a -= a_offset;

    --work;



    /* Function Body */

    ifault = 1;



    i = n;

    for (i1 = 2; i1 <= n; ++i1) {

	l  = i-2;

	f  = a[i + (i-1)*ndima];

	gl = 0;



	if (l < 1) goto L25;



	for (k = 1; k <= l; ++k) {

	    d__1 = a[i + k*ndima];

	    gl  += d__1*d__1;

	}

L25:

	h = gl + f*f;



	if (gl > 1e-35) goto L30;



	work[i]     = 0;

	work[n + i] = f;

	goto L65;

L30:

	++l;

	gl = TMath_Sqrt(h);

	if (f >= 0) gl = -gl;

	work[n + i] = gl;

	h -= f*gl;

	a[i + (i-1)*ndima] = f - gl;

	f = 0;

	for (j = 1; j <= l; ++j) {

	    a[j + i*ndima] = a[i + j*ndima] / h;

	    gl = 0;

	    for (k = 1; k <= j; ++k) { gl += a[j + k*ndima]*a[i + k*ndima]; }

	    if (j >= l) goto L47;

	    j1 = j + 1;

	    for (k = j1; k <= l; ++k) {	gl += a[k + j*ndima]*a[i + k*ndima]; }

L47:

	    work[n + j] = gl / h;

	    f += gl*a[j + i*ndima];

	}

	hh = f / (h + h);

	for (j = 1; j <= l; ++j) {

	    f  = a[i + j*ndima];

	    gl = work[n + j] - hh*f;

	    work[n + j] = gl;

	    for (k = 1; k <= j; ++k) {

		a[j + k*ndima] = a[j + k*ndima] - f*work[n + k] - gl*a[i + k*ndima];

	    }

	}

	work[i] = h;

L65:

	--i;

    }

    work[1] = 0;

    work[n + 1] = 0;

    for (i = 1; i <= n; ++i) {

	l = i-1;

	if (work[i] == 0 || l == 0) goto L100;



	for (j = 1; j <= l; ++j) {

	    gl = 0;

	    for (k = 1; k <= l; ++k) { gl += a[i + k*ndima]*a[k + j*ndima]; }

	    for (k = 1; k <= l; ++k) { a[k + j*ndima] -= gl*a[k + i*ndima]; }

	}

L100:

	work[i] = a[i + i*ndima];

	a[i + i*ndima] = 1;

	if (l == 0) continue;



	for (j = 1; j <= l; ++j) {

	    a[i + j*ndima] = 0;

	    a[j + i*ndima] = 0;

	}

    }



    n1 = n - 1;

    for (i = 2; i <= n; ++i) {

	i0 = n + i-1;

	work[i0] = work[i0 + 1];

    }

    work[n + n] = 0;

    b = 0;

    f = 0;

    for (l = 1; l <= n; ++l) {

	j = 0;

	h = precis*(TMath_Abs(work[l]) + TMath_Abs(work[n + l]));

	if (b < h) b = h;

	for (m1 = l; m1 <= n; ++m1) {

	    m = m1;

	    if (TMath_Abs(work[n + m]) <= b)	goto L150;

	}



L150:

	if (m == l) goto L205;



L160:

	if (j == mits) return;

	++j;

	pt = (work[l + 1] - work[l]) / (work[n + l]*2);

	r  = TMath_Sqrt(pt*pt + 1);

	pr = pt + r;

	if (pt < 0) pr = pt - r;



	h = work[l] - work[n + l] / pr;

	for (i = l; i <= n; ++i) { work[i] -= h; }

	f += h;

	pt = work[m];

	c  = 1;

	s  = 0;

	m1 = m - 1;

	i  = m;

	for (i1 = l; i1 <= m1; ++i1) {

	    j = i;

	    --i;

	    gl = c*work[n + i];

	    h  = c*pt;

	    if (TMath_Abs(pt) >= TMath_Abs(work[n + i])) goto L180;



	    c = pt / work[n + i];

	    r = TMath_Sqrt(c*c + 1);

	    work[n + j] = s*work[n + i]*r;

	    s  = 1 / r;

	    c /= r;

	    goto L190;

L180:

	    c = work[n + i] / pt;

	    r = TMath_Sqrt(c*c + 1);

	    work[n + j] = s*pt*r;

	    s = c / r;

	    c = 1 / r;

L190:

	    pt = c*work[i] - s*gl;

	    work[j] = h + s*(c*gl + s*work[i]);

	    for (k = 1; k <= n; ++k) {

		h = a[k + j*ndima];

		a[k + j*ndima] = s*a[k + i*ndima] + c*h;

		a[k + i*ndima] = c*a[k + i*ndima] - s*h;

	    }

	}

	work[n + l] = s*pt;

	work[l]     = c*pt;



	if (TMath_Abs(work[n + l]) > b) goto L160;



L205:

	work[l] += f;

    }

    for (i = 1; i <= n1; ++i) {

	k  = i;

	pt = work[i];

	i1 = i + 1;

	for (j = i1; j <= n; ++j) {

	    if (work[j] >= pt) continue;

	    k  = j;

	    pt = work[j];

	}



	if (k == i) continue;



	work[k] = work[i];

	work[i] = pt;

	for (j = 1; j <= n; ++j) {

	    pt = a[j + i*ndima];

	    a[j + i*ndima] = a[j + k*ndima];

	    a[j + k*ndima] = pt;

	}

    }

    ifault = 0;

} /* mneig_ */



//______________________________________________________________________________

void Midnight::mnemat(MDouble *emat, MInt ndim)

{

//*-*-*-*-*-*Calculates the external error matrix from the internal matrix*-*

//*-*        =============================================================

//*-*        Calculates the external error matrix from the internal

//*-*        to be called by user.

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    /* System generated locals */

    MInt emat_dim1, emat_offset;



    /* Local variables */

    static MDouble dxdi, dxdj;

    static MInt i, j, k, npard, k2, kk, iz, nperln, kga, kgb;

    static MString ctemp;



    /* Parameter adjustments */

    emat_dim1 = ndim;

    emat_offset = emat_dim1 + 1;

    emat -= emat_offset;



    /* Function Body */

    if (fISW[1] < 1) return;

    if (fISW[4] >= 2) {

	fPrintf(" EXTERNAL ERROR MATRIX.    NDIM=%4d    NPAR=%3d    ERR DEF=%g",ndim,fNpar,fUp);

    }

//*-*-                   size of matrix to be printed

    npard = fNpar;

    if (ndim < fNpar) {

	npard = ndim;

	if (fISW[4] >= 0) {

	    fPrintf(" USER-DIMENSIONED  ARRAY EMAT NOT BIG ENOUGH. REDUCED MATRIX CALCULATED.");

	}

    }

//*-*-                NPERLN is the number of elements that fit on one line



    nperln = (fNpagwd - 5) / 10;

    nperln = TMath_Min(nperln,13);

    if (fISW[4] >= 1 && npard > nperln) {

	fPrintf(" ELEMENTS ABOVE DIAGONAL ARE NOT PRINTED.");

    }

//*-*-                I counts the rows of the matrix

    for (i = 1; i <= npard; ++i) {

	mndxdi(fX[i-1], i-1, dxdi);

	kga = i*(i-1) / 2;

	for (j = 1; j <= i; ++j) {

	    mndxdi(fX[j-1], j-1, dxdj);

	    kgb = kga + j;

	    emat[i + j*emat_dim1] = dxdi*fVhmat[kgb-1]*dxdj*fUp;

	    emat[j + i*emat_dim1] = emat[i + j*emat_dim1];

	}

    }

//*-*-                   IZ is number of columns to be printed in row I

    if (fISW[4] >= 2) {

	for (i = 1; i <= npard; ++i) {

	    iz = npard;

	    if (npard >= nperln) iz = i;

            ctemp = " ";

	    for (k = 1; nperln < 0 ? k >= iz : k <= iz; k += nperln) {

		k2 = k + nperln - 1;

		if (k2 > iz) k2 = iz;

		for (kk = k; kk <= k2; ++kk) {

		    ctemp += Form("%10.3e ",emat[i + kk*emat_dim1]);

	        }

                fPrintf("%s",(const char*)ctemp);

	    }

	}

    }

} /* mnemat_ */



//______________________________________________________________________________

void Midnight::mnerrs(MInt number, MDouble &eplus, MDouble &eminus, MDouble &eparab, MDouble &gcc)

{

//*-*-*-*-*-*-*-*-*-*Utility routine to get MINOS errors*-*-*-*-*-*-*-*-*-*-*

//*-*                ===================================

//*-*    Called by user.

//*-*    NUMBER is the parameter number

//*-*    values returned by MNERRS:

//*-*       EPLUS, EMINUS are MINOS errors of parameter NUMBER,

//*-*       EPARAB is 'parabolic' error (from error matrix).

//*-*                 (Errors not calculated are set = 0)

//*-*       GCC is global correlation coefficient from error matrix

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    static MDouble dxdi;

    static MInt ndiag, iin, iex;



    iex = number+1;



    if (iex > fNu || iex <= 0) goto L900;

    iin = fNiofex[iex-1];

    if (iin <= 0) goto L900;



//*-*-            IEX is external number, IIN is internal number

    eplus  = fErp[iin-1];

    if (eplus == fUndefi)  eplus = 0;

    eminus = fErn[iin-1];

    if (eminus == fUndefi) eminus = 0;

    mndxdi(fX[iin-1], iin-1, dxdi);

    ndiag  = iin*(iin + 1) / 2;

    eparab = TMath_Abs(dxdi*TMath_Sqrt(TMath_Abs(fUp*fVhmat[ndiag- 1])));

//*-*-             global correlation coefficient

    gcc = 0;

    if (fISW[1] < 2) return;

    gcc = fGlobcc[iin-1];

    return;

//*-*-                 ERROR.  parameter number not valid

L900:

    eplus  = 0;

    eminus = 0;

    eparab = 0;

    gcc    = 0;

} /* mnerrs_ */



//______________________________________________________________________________

void Midnight::mneval(MDouble anext, MDouble &fnext, MInt &ierev)

{

//*-*-*-*-*-*-*Evaluates the function being analyzed by MNCROS*-*-*-*-*-*-*-*

//*-*          ===============================================

//*-*      Evaluates the function being analyzed by MNCROS, which is

//*-*      generally the minimum of FCN with respect to all remaining

//*-*      variable parameters.  The class data members contains the

//*-*      data necessary to know the values of U(KE1CR) and U(KE2CR)

//*-*      to be used, namely     U(KE1CR) = XMIDCR + ANEXT*XDIRCR

//*-*      and (if KE2CR .NE. 0)  U(KE2CR) = YMIDCR + ANEXT*YDIRCR

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    static MInt nparx;



    fU[fKe1cr-1] = fXmidcr + anext*fXdircr;

    if (fKe2cr != 0) fU[fKe2cr-1] = fYmidcr + anext*fYdircr;

    mninex(fX);

    nparx = fNpar;

    (*fFCN)(nparx, fGin, fnext, fU, 4);    ++fNfcn;

    ierev = 0;

    if (fNpar > 0) {

	fItaur = 1;

	fAmin = fnext;

	fISW[0] = 0;

	mnmigr();

	fItaur = 0;

	fnext = fAmin;

	if (fISW[0] >= 1) ierev = 1;

	if (fISW[3] < 1)  ierev = 2;

    }

} /* mneval_ */



//______________________________________________________________________________

void Midnight::mnexcm(MString comand, MDouble *plist, MInt llist, MInt &ierflg)

{

//*-*-*-*-*-*Interprets a command and takes appropriate action*-*-*-*-*-*-*-*

//*-*        =================================================

//*-*        either directly by skipping to the corresponding code in

//*-*        MNEXCM, or by setting up a call to a function

//*-*

//*-*  recognized MINUIT commands:

//*-*  obsolete commands:

//*-*      IERFLG is now (94.5) defined the same as ICONDN in MNCOMD

//*-*            = 0: command executed normally

//*-*              1: command is blank, ignored

//*-*              2: command line unreadable, ignored

//*-*              3: unknown command, ignored

//*-*              4: abnormal termination (e.g., MIGRAD not converged)

//*-*              9: reserved

//*-*             10: END command

//*-*             11: EXIT or STOP command

//*-*             12: RETURN command

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    /* Initialized data */



    static MString clower = "abcdefghijklmnopqrstuvwxyz";

    static MString cupper = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

    static MString cname[40] = {

      "MINImize  ",

      "SEEk      ",

      "SIMplex   ",

      "MIGrad    ",

      "MINOs     ",

      "SET xxx   ",

      "SHOw xxx  ",

      "TOP of pag",

      "FIX       ",

      "REStore   ",

      "RELease   ",

      "SCAn      ",

      "CONtour   ",

      "HESse     ",

      "SAVe      ",

      "IMProve   ",

      "CALl fcn  ",

      "STAndard  ",

      "END       ",

      "EXIt      ",

      "RETurn    ",

      "CLEar     ",

      "HELP      ",

      "MNContour ",

      "STOp      ",

      "JUMp      ",

      "          ",

      "          ",

      "          ",

      "          ",

      "          ",

      "          ",

      "          ",

      "COVARIANCE",

      "PRINTOUT  ",

      "GRADIENT  ",

      "MATOUT    ",

      "ERROR DEF ",

      "LIMITS    ",

      "PUNCH     "};



    static MInt nntot = 40;



    /* Local variables */

    static MDouble step, xptu[101], yptu[101], f, rno;

    static MInt icol, kcol, ierr, iint, iext, lnow, nptu, i, iflag, ierrf;

    static MInt ilist, nparx, izero, nf, lk, it, iw, inonde, nsuper;

    static MInt it2, ke1, ke2, nowprt, kll, let, krl;

    static MString chwhy, c26, cvblnk, cneway, comd;

    static MString ctemp;

    static MBool lfreed, ltofix, lfixed;



//*-*  alphabetical order of command names!



    /* Function Body */



    lk = strlen((const char*)comand);

    if (lk > 20) lk = 20;

    fCword =  comand;

//*-*-             get upper case

    for (icol = 1; icol <= lk; ++icol) {

	for (let = 1; let <= 26; ++let) {

	    if (fCword[icol-1] == clower[let-1]) {

		fCword[icol-1] = cupper[let-1];

	    }

	}

    }

//*-*-          Copy the first MAXP arguments into WORD7, making

//*-*-          sure that WORD7(1)=0 if LLIST=0

    for (iw = 1; iw <= 30; ++iw) {

	fWord7[iw-1] = 0;

	if (iw <= llist) fWord7[iw-1] = plist[iw-1];

    }

    ++fIcomnd;

    fNfcnlc = fNfcn;

    if (fCword(0,7) != "SET PRI" || fWord7[0] >= 0) {

	if (fISW[4] >= 0) {

	    lnow = llist;

	    if (lnow > 4) lnow = 4;

	    fPrintf(" **********");

            ctemp = Form(" **%5d **%s",fIcomnd,(const char*)fCword);

	    for (i = 1; i <= lnow; ++i) {

                ctemp += Form("%12.4g",plist[i-1]);

	    }

            fPrintf("%s",(const char*)ctemp);

	    inonde = 0;

	    if (llist > lnow) {

		kll = llist;

		if (llist > 30) {

		    inonde = 1;

		    kll = 30;

		}

                fPrintf(" ***********");

		for (i = lnow + 1; i <= kll; ++i) {

                   fPrintf("%12.4g",plist[i-1]);

		}

	    }

	    fPrintf(" **********");

	    if (inonde > 0) {

		fPrintf("  ERROR: ABOVE CALL TO MNEXCM TRIED TO PASS MORE THAN 30 PARAMETERS.");

	    }

	}

    }

    fNfcnmx = MInt(fWord7[0]);

    if (fNfcnmx <= 0) {

	fNfcnmx = fNpar*100 + 200 + fNpar*fNpar*5;

    }

    fEpsi = fWord7[1];

    if (fEpsi <= 0) {

	fEpsi = fUp*.1;

    }

    fLnewmn = kFALSE;

    fLphead = kTRUE;

    fISW[0] = 0;

    ierflg = 0;

//*-*-               look for command in list CNAME . . . . . . . . . .

    for (i = 1; i <= nntot; ++i) {

	if (fCword(0,3) == cname[i-1](0,3)) goto L90;

    }

    fPrintf("UNKNOWN COMMAND IGNORED:%s", (const char*)comand);

    ierflg = 3;

    return;

//*-*-               normal case: recognized MINUIT command . . . . . . .

L90:

    if (fCword(0,4) == "MINO") i = 5;

    if (i != 6 && i != 7 && i != 8 && i != 23) {

	fCfrom  = cname[i-1];

	fNfcnfr = fNfcn;

    }

//*-*-             1    2    3    4    5    6    7    8    9   10

    switch (i) {

	case 1:  goto L400;

	case 2:  goto L200;

	case 3:  goto L300;

	case 4:  goto L400;

	case 5:  goto L500;

	case 6:  goto L700;

	case 7:  goto L700;

	case 8:  goto L800;

	case 9:  goto L900;

	case 10:  goto L1000;

	case 11:  goto L1100;

	case 12:  goto L1200;

	case 13:  goto L1300;

	case 14:  goto L1400;

	case 15:  goto L1500;

	case 16:  goto L1600;

	case 17:  goto L1700;

	case 18:  goto L1800;

	case 19:  goto L1900;

	case 20:  goto L1900;

	case 21:  goto L1900;

	case 22:  goto L2200;

	case 23:  goto L2300;

	case 24:  goto L2400;

	case 25:  goto L1900;

	case 26:  goto L2600;

	case 27:  goto L3300;

	case 28:  goto L3300;

	case 29:  goto L3300;

	case 30:  goto L3300;

	case 31:  goto L3300;

	case 32:  goto L3300;

	case 33:  goto L3300;

	case 34:  goto L3400;

	case 35:  goto L3500;

	case 36:  goto L3600;

	case 37:  goto L3700;

	case 38:  goto L3800;

	case 39:  goto L3900;

	case 40:  goto L4000;

    }

//*-*-                                       . . . . . . . . . . seek

L200:

    mnseek();

    return;

//*-*-                                       . . . . . . . . . . simplex

L300:

    mnsimp();

    if (fISW[3] < 1) ierflg = 4;

    return;

//*-*-                                       . . . . . . migrad, minimize

L400:

    nf = fNfcn;

    fApsi = fEpsi;

    mnmigr();

    mnwerr();

    if (fISW[3] >= 1) return;

    ierflg = 4;

    if (fISW[0] == 1) return;

    if (fCword(0,3) == "MIG") return;



    fNfcnmx = fNfcnmx + nf - fNfcn;

    nf = fNfcn;

    mnsimp();

    if (fISW[0] == 1) return;

    fNfcnmx = fNfcnmx + nf - fNfcn;

    mnmigr();

    if (fISW[3] >= 1) ierflg = 0;

    mnwerr();

    return;

//*-*-                                       . . . . . . . . . . minos

L500:

    nsuper = fNfcn + ((fNpar + 1) << 1)*fNfcnmx;

//*-*-         possible loop over new minima

    fEpsi = fUp*.1;

L510:

    mncuve();

    mnmnos();

    if (! fLnewmn) return;

    mnrset(0);

    mnmigr();

    mnwerr();

    if (fNfcn < nsuper) goto L510;

    fPrintf(" TOO MANY FUNCTION CALLS. MINOS GIVES UP");

    ierflg = 4;

    return;

//*-*-                                       . . . . . . . . . .set, show

L700:

    mnset();

    return;

//*-*-                                       . . . . . . . . . . top of page



L800:

    fPrintf("1");

    return;

//*-*-                                       . . . . . . . . . . fix

L900:

    ltofix = kTRUE;

//*-*-                                       . . (also release) ....

L901:

    lfreed = kFALSE;

    lfixed = kFALSE;

    if (llist == 0) {

	fPrintf("%s:  NO PARAMETERS REQUESTED ",(const char*)fCword);

	return;

    }

    for (ilist = 1; ilist <= llist; ++ilist) {

	iext = MInt(plist[ilist-1]);

	chwhy = " IS UNDEFINED.";

	if (iext <= 0) goto L930;

	if (iext > fNu) goto L930;

	if (fNvarl[iext-1] < 0) goto L930;

	chwhy = " IS CONSTANT.  ";

	if (fNvarl[iext-1] == 0) goto L930;

	iint = fNiofex[iext-1];

	if (ltofix) {

	    chwhy = " ALREADY FIXED.";

	    if (iint == 0) goto L930;

	    mnfixp(iint-1, ierr);

	    if (ierr == 0) lfixed = kTRUE;

	    else           ierflg = 4;

	} else {

	    chwhy = " ALREADY VARIABLE.";

	    if (iint > 0) goto L930;

	    krl = -abs(iext);

	    mnfree(krl);

	    lfreed = kTRUE;

	}

	continue;

L930:

	fPrintf(" PARAMETER%4d %s IGNORED.",iext,(const char*)chwhy);

    }

    if (lfreed || lfixed) mnrset(0);

    if (lfreed) {

	fISW[1] = 0;

	fDcovar = 1;

	fEDM = fBigedm;

	fISW[3] = 0;

    }

    mnwerr();

    if (fISW[4] > 1) mnprin(5, fAmin);

    return;

//*-*-                                       . . . . . . . . . . restore

L1000:

    it = MInt(fWord7[0]);

    if (it > 1 || it < 0) goto L1005;

    lfreed = fNpfix > 0;

    mnfree(it);

    if (lfreed) {

	mnrset(0);

	fISW[1] = 0;

	fDcovar = 1;

	fEDM    = fBigedm;

    }

    return;

L1005:

    fPrintf(" IGNORED.  UNKNOWN ARGUMENT:%4d",it);

    ierflg = 3;

    return;

//*-*-                                       . . . . . . . . . . release

L1100:

    ltofix = kFALSE;

    goto L901;

//*-*-                                      . . . . . . . . . . scan . . .

L1200:

    iext = MInt(fWord7[0]);

    if (iext <= 0) goto L1210;

    it2 = 0;

    if (iext <= fNu) it2 = fNiofex[iext-1];

    if (it2 <= 0) goto L1250;



L1210:

    mnscan();

    return;

L1250:

    fPrintf(" PARAMETER%4d NOT VARIABLE.",iext);

    ierflg = 3;

    return;

//*-*-                                       . . . . . . . . . . contour

L1300:

    ke1 = MInt(fWord7[0]);

    ke2 = MInt(fWord7[1]);

    if (ke1 == 0) {

	if (fNpar == 2) {

	    ke1 = fNexofi[0];

	    ke2 = fNexofi[1];

	} else {

	    fPrintf("%s:  NO PARAMETERS REQUESTED ",(const char*)fCword);

	    ierflg = 3;

	    return;

	}

    }

    fNfcnmx = 1000;

    mncntr(ke1, ke2, ierrf);

    if (ierrf > 0) ierflg = 3;

    return;

//*-*-                                       . . . . . . . . . . hesse

L1400:

    mnhess();

    mnwerr();

    if (fISW[4] >= 0) mnprin(2, fAmin);

    if (fISW[4] >= 1) mnmatu(1);

    return;

//*-*-                                       . . . . . . . . . . save

L1500:

    mnsave();

    return;

//*-*-                                       . . . . . . . . . . improve

L1600:

    mncuve();

    mnimpr();

    if (fLnewmn) goto L400;

    ierflg = 4;

    return;

//*-*-                                       . . . . . . . . . . call fcn

L1700:

    iflag = MInt(fWord7[0]);

    nparx = fNpar;

    f = fUndefi;

    (*fFCN)(nparx, fGin, f, fU, iflag);    ++fNfcn;

    nowprt = 0;

    if (f != fUndefi) {

	if (fAmin == fUndefi) {

	    fAmin  = f;

	    nowprt = 1;

	} else if (f < fAmin) {

	    fAmin  = f;

	    nowprt = 1;

	}

	if (fISW[4] >= 0 && iflag <= 5 && nowprt == 1) {

	    mnprin(5, fAmin);

	}

	if (iflag == 3)  fFval3 = f;

    }

    if (iflag > 5) mnrset(1);

    return;

//*-*-                                       . . . . . . . . . . standard

L1800:

//    stand();

    return;

//*-*-                                      . . . return, stop, end, exit

L1900:

    it = MInt(fWord7[0]);

    if (fFval3 != fAmin && it == 0) {

	iflag = 3;

	fPrintf(" CALL TO USER FUNCTION WITH IFLAG = 3");

	nparx = fNpar;

	(*fFCN)(nparx, fGin, f, fU, iflag);	++fNfcn;

    }

    ierflg = 11;

    if (fCword(0,3) == "END") ierflg = 10;

    if (fCword(0,3) == "RET") ierflg = 12;

    return;

//*-*-                                       . . . . . . . . . . clear

L2200:

    mncler();

    if (fISW[4] >= 1) {

	fPrintf(" MINUIT MEMORY CLEARED. NO PARAMETERS NOW DEFINED.");

    }

    return;

//*-*-                                       . . . . . . . . . . help

L2300:

    kcol = 0;

    for (icol = 5; icol <= lk; ++icol) {

	if (fCword[icol-1] == ' ') continue;

	kcol = icol;

	goto L2320;

    }

L2320:

    if (kcol == 0) comd = "*   ";

    else           comd = fCword[kcol-1];

    mnhelp(comd);

    return;

//*-*-                                      . . . . . . . . . . MNContour

L2400:

    fEpsi = fUp*.05;

    ke1 = MInt(fWord7[0]);

    ke2 = MInt(fWord7[1]);

    if (ke1 == 0 && fNpar == 2) {

	ke1 = fNexofi[0];

	ke2 = fNexofi[1];

    }

    nptu = MInt(fWord7[2]);

    if (nptu <= 0)  nptu = 20;

    if (nptu > 101) nptu = 101;

    fNfcnmx = (nptu + 5)*100*(fNpar + 1);

    mncont(ke1, ke2, nptu, xptu, yptu, ierrf);

    if (ierrf < nptu) ierflg = 4;

    if (ierrf == -1)  ierflg = 3;

    return;

//*-*-                                     . . . . . . . . . . jump

L2600:

    step = fWord7[0];

    if (step <= 0) step = 2;

    rno = 0;

    izero = 0;

    for (i = 1; i <= fNpar; ++i) {

	mnrn15(rno, izero);

	rno      = rno*2 - 1;

	fX[i-1] += rno*step*fWerr[i-1];

    }

    mninex(fX);

    mnamin();

    mnrset(0);

    return;

//*-*-                                     . . . . . . . . . . blank line

L3300:

    fPrintf(" BLANK COMMAND IGNORED.");

    ierflg = 1;

    return;

//*-*  . . . . . . . . obsolete commands     . . . . . . . . . . . . . .

//*-*-                                     . . . . . . . . . . covariance

L3400:

    fPrintf(" THE *COVARIANCE* COMMAND IS OSBSOLETE. THE COVARIANCE MATRIX IS NOW SAVED IN A DIFFERENT FORMAT WITH THE *SAVE* COMMAND AND READ IN WITH:*SET COVARIANCE*");

    ierflg = 3;

    return;

//*-*-                                       . . . . . . . . . . printout

L3500:

    cneway = "SET PRInt ";

    goto L3100;

//*-*-                                       . . . . . . . . . . gradient

L3600:

    cneway = "SET GRAd  ";

    goto L3100;

//*-*-                                       . . . . . . . . . . matout

L3700:

    cneway = "SHOW COVar";

    goto L3100;

//*-*-                                       . . . . . . . . . error def

L3800:

    cneway = "SET ERRdef";

    goto L3100;

//*-*-                                       . . . . . . . . . . limits

L3900:

    cneway = "SET LIMits";

    goto L3100;

//*-*-                                       . . . . . . . . . . punch

L4000:

    cneway = "SAVE      ";

//*-*-                               ....... come from obsolete commands

L3100:

    fPrintf(" OBSOLETE COMMAND:%s   PLEASE USE: %s",(const char*)fCword

                                                 ,(const char*)cneway);

    fCword = cneway;

    if (fCword == "SAVE      ") goto L1500;

    goto L700;

//*-*                                 . . . . . . . . . . . . . . . . . .

} /* mnexcm_ */



//______________________________________________________________________________

void Midnight::mnexin(MDouble *pint)

{

//*-*-*-*-*Transforms the external parameter values U to internal values*-*-*

//*-*      =============================================================

//*-*        Transforms the external parameter values U to internal

//*-*        values in the dense array PINT.

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    static MDouble pinti;

    static MInt iint, iext;



    fLimset = kFALSE;

    for (iint = 1; iint <= fNpar; ++iint) {

	iext = fNexofi[iint-1];

	mnpint(fU[iext-1], iext-1, pinti);

	pint[iint-1] = pinti;

    }

} /* mnexin_ */



//______________________________________________________________________________

void Midnight::mnfixp(MInt iint1, MInt &ierr)

{

//*-*-*-*-*-*-*Removes parameter IINT from the internal parameter list*-*-*

//*-*          =======================================================

//*-*        and arranges the rest of the list to fill the hole.

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    /* Local variables */

    static MDouble yy[kMAXDIM], yyover;

    static MInt kold, nold, ndex, knew, iext, i, j, m, n, lc, ik;



//*-*-                          first see if it can be done

    ierr = 0;

    MInt iint = iint1+1;

    if (iint > fNpar || iint <= 0) {

	ierr = 1;

        fPrintf(" MINUIT ERROR.  ARGUMENT TO MNFIXP=%4d",iint);

	return;

    }

    iext = fNexofi[iint-1];

    if (fNpfix >= fMaxpar) {

	ierr = 1;

	fPrintf(" MINUIT CANNOT FIX PARAMETER%4d MAXIMUM NUMBER THAT CAN BE FIXED IS %d",iext,fMaxpar);

	return;

    }

//*-*-                          reduce number of variable parameters by one



    fNiofex[iext-1] = 0;

    nold = fNpar;

    --fNpar;

//*-*-                      save values in case parameter is later restored



    ++fNpfix;

    fIpfix[fNpfix-1]  = iext;

    lc                = iint;

    fXs[fNpfix-1]     = fX[lc-1];

    fXts[fNpfix-1]    = fXt[lc-1];

    fDirins[fNpfix-1] = fWerr[lc-1];

    fGrds[fNpfix-1]   = fGrd[lc-1];

    fG2s[fNpfix-1]    = fG2[lc-1];

    fGsteps[fNpfix-1] = fGstep[lc-1];

//*-*-                       shift values for other parameters to fill hole

    for (ik = iext + 1; ik <= fNu; ++ik) {

	if (fNiofex[ik-1] > 0) {

	    lc = fNiofex[ik-1] - 1;

	    fNiofex[ik-1] = lc;

	    fNexofi[lc-1] = ik;

	    fX[lc-1]      = fX[lc];

	    fXt[lc-1]     = fXt[lc];

	    fDirin[lc-1]  = fDirin[lc];

	    fWerr[lc-1]   = fWerr[lc];

	    fGrd[lc-1]    = fGrd[lc];

	    fG2[lc-1]     = fG2[lc];

	    fGstep[lc-1]  = fGstep[lc];

	}

    }

    if (fISW[1] <= 0) return;

//*-*-                   remove one row and one column from variance matrix

    if (fNpar <= 0)   return;

    for (i = 1; i <= nold; ++i) {

	m       = TMath_Max(i,iint);

	n       = TMath_Min(i,iint);

	ndex    = m*(m-1) / 2 + n;

	yy[i-1] = fVhmat[ndex-1];

    }

    yyover = 1 / yy[iint-1];

    knew   = 0;

    kold   = 0;

    for (i = 1; i <= nold; ++i) {

	for (j = 1; j <= i; ++j) {

	    ++kold;

	    if (j == iint || i == iint) continue;

	    ++knew;

	    fVhmat[knew-1] = fVhmat[kold-1] - yy[j-1]*yy[i-1]*yyover;

	}

    }

} /* mnfixp_ */



//______________________________________________________________________________

void Midnight::mnfree(MInt k)

{

//*-*-*-*Restores one or more fixed parameter(s) to variable status*-*-*-*-*-*

//*-*    ==========================================================

//*-*        Restores one or more fixed parameter(s) to variable status

//*-*        by inserting it into the internal parameter list at the

//*-*        appropriate place.

//*-*

//*-*        K = 0 means restore all parameters

//*-*        K = 1 means restore the last parameter fixed

//*-*        K = -I means restore external parameter I (if possible)

//*-*        IQ = fix-location where internal parameters were stored

//*-*        IR = external number of parameter being restored

//*-*        IS = internal number of parameter being restored

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    /* Local variables */

    static MDouble grdv, xv, dirinv, g2v, gstepv, xtv;

    static MInt i, ipsav, ka, lc, ik, iq, ir, is;



    if (k > 1) {

	fPrintf(" CALL TO MNFREE IGNORED.  ARGUMENT GREATER THAN ONE");

    }

    if (fNpfix < 1) {

	fPrintf(" CALL TO MNFREE IGNORED.  THERE ARE NO FIXED PARAMETERS");

    }

    if (k == 1 || k == 0) goto L40;



//*-*-                  release parameter with specified external number

    ka = abs(k);

    if (fNiofex[ka-1] == 0) goto L15;

    fPrintf(" IGNORED.  PARAMETER SPECIFIED IS ALREADY VARIABLE.");

    return;

L15:

    if (fNpfix < 1) goto L21;

    for (ik = 1; ik <= fNpfix; ++ik) { if (fIpfix[ik-1] == ka) goto L24; }

L21:

    fPrintf(" PARAMETER%4d NOT FIXED.  CANNOT BE RELEASED.",ka);

    return;

L24:

    if (ik == fNpfix) goto L40;



//*-*-                  move specified parameter to end of list

    ipsav  = ka;

    xv     = fXs[ik-1];

    xtv    = fXts[ik-1];

    dirinv = fDirins[ik-1];

    grdv   = fGrds[ik-1];

    g2v    = fG2s[ik-1];

    gstepv = fGsteps[ik-1];

    for (i = ik + 1; i <= fNpfix; ++i) {

	fIpfix[i-2]  = fIpfix[i-1];

	fXs[i-2]     = fXs[i-1];

	fXts[i-2]    = fXts[i-1];

	fDirins[i-2] = fDirins[i-1];

	fGrds[i-2]   = fGrds[i-1];

	fG2s[i-2]    = fG2s[i-1];

	fGsteps[i-2] = fGsteps[i-1];

    }

    fIpfix[fNpfix-1]  = ipsav;

    fXs[fNpfix-1]     = xv;

    fXts[fNpfix-1]    = xtv;

    fDirins[fNpfix-1] = dirinv;

    fGrds[fNpfix-1]   = grdv;

    fG2s[fNpfix-1]    = g2v;

    fGsteps[fNpfix-1] = gstepv;

//*-*-               restore last parameter in fixed list  -- IPFIX(NPFIX)

L40:

    if (fNpfix < 1) goto L300;

    ir = fIpfix[fNpfix-1];

    is = 0;

    for (ik = fNu; ik >= ir; --ik) {

	if (fNiofex[ik-1] > 0) {

	    lc = fNiofex[ik-1] + 1;

	    is = lc - 1;

	    fNiofex[ik-1] = lc;

	    fNexofi[lc-1] = ik;

	    fX[lc-1]      = fX[lc-2];

	    fXt[lc-1]     = fXt[lc-2];

	    fDirin[lc-1]  = fDirin[lc-2];

	    fWerr[lc-1]   = fWerr[lc-2];

	    fGrd[lc-1]    = fGrd[lc-2];

	    fG2[lc-1]     = fG2[lc-2];

	    fGstep[lc-1]  = fGstep[lc-2];

	}

    }

    ++fNpar;

    if (is == 0) is = fNpar;

    fNiofex[ir-1] = is;

    fNexofi[is-1] = ir;

    iq           = fNpfix;

    fX[is-1]     = fXs[iq-1];

    fXt[is-1]    = fXts[iq-1];

    fDirin[is-1] = fDirins[iq-1];

    fWerr[is-1]  = fDirins[iq-1];

    fGrd[is-1]   = fGrds[iq-1];

    fG2[is-1]    = fG2s[iq-1];

    fGstep[is-1] = fGsteps[iq-1];

    --fNpfix;

    fISW[1] = 0;

    fDcovar = 1;

    if (fISW[4] - fItaur >= 1) {

	fPrintf("                    PARAMETER%4d  %s RESTORED TO VARIABLE.",ir,

                      (const char*)fCpnam[ir-1]);

    }

    if (k == 0) goto L40;

L300:

//*-*-        if different from internal, external values are taken

    mnexin(fX);

} /* mnfree_ */



//______________________________________________________________________________

void Midnight::mngrad()

{

//*-*-*-*-*-*-*-*-*-*Interprets the SET GRAD command*-*-*-*-*-*-*-*-*-*-*-*-*

//*-*                ===============================

//*-*       Called from MNSET

//*-*       Interprets the SET GRAD command, which informs MINUIT whether

//*-*       the first derivatives of FCN will be calculated by the user

//*-*       inside FCN.  It can check the user derivative calculation

//*-*       by comparing it with a finite difference approximation.

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    /* Local variables */

    static MDouble gf[kMAXDIM], fzero, err;

    static MInt i, nparx, lc, istsav;

    static MBool lnone;

    static MString cwd = "    ";



    fISW[2] = 1;

    nparx   = fNpar;

    if (fWord7[0] > 0) goto L2000;



//*-*-                 get user-calculated first derivatives from FCN

    for (i = 1; i <= fNu; ++i) { fGin[i-1] = fUndefi; }

    mninex(fX);

    (*fFCN)(nparx, fGin, fzero, fU, 2);    ++fNfcn;

    mnderi();

    for (i = 1; i <= fNpar; ++i) { gf[i-1] = fGrd[i-1]; }

//*-*-                   get MINUIT-calculated first derivatives

    fISW[2] = 0;

    istsav  = fIstrat;

    fIstrat = 2;

    mnhes1();

    fIstrat = istsav;

    fPrintf(" CHECK OF GRADIENT CALCULATION IN FCN");

    fPrintf("            PARAMETER      G(IN FCN)   G(MINUIT)  DG(MINUIT)   AGREEMENT");

    fISW[2] = 1;

    lnone = kFALSE;

    for (lc = 1; lc <= fNpar; ++lc) {

	i   = fNexofi[lc-1];

	cwd = "GOOD";

	err = fDgrd[lc-1];

	if (TMath_Abs(gf[lc-1] - fGrd[lc-1]) > err)  cwd = " BAD";

	if (fGin[i-1] == fUndefi) {

	    cwd      = "NONE";

	    lnone    = kTRUE;

	    gf[lc-1] = 0;

	}

	if (cwd != "GOOD") fISW[2] = 0;

	fPrintf("       %5d  %10s%12.4e%12.4e%12.4e    %s",i

                    ,(const char*)fCpnam[i-1]

                    ,gf[lc-1],fGrd[lc-1],err,(const char*)cwd);

    }

    if (lnone) {

	fPrintf("  AGREEMENT=NONE  MEANS FCN DID NOT CALCULATE THE DERIVATIVE");

    }

    if (fISW[2] == 0) {

	fPrintf(" MINUIT DOES NOT ACCEPT DERIVATIVE CALCULATIONS BY FCN");

        fPrintf(" TO FORCE ACCEPTANCE, ENTER *SET GRAD    1*");

    }



L2000:

    return;

} /* mngrad_ */



//______________________________________________________________________________

void Midnight::mnhelp(MString comd)

{

//*-*-*-*-*-*-*-*HELP routine for MINUIT interactive commands*-*-*-*-*-*-*-*-*

//*-*            ============================================

//*-*

//*-*      COMD ='*   '  prints a global help for all commands

//*-*      COMD =Command_name: print detailed help for one command.

//*-*         Note that at least 3 characters must be given for the command

//*-*         name.

//*-*

//*-*     Author: Rene Brun

//*-*             comments extracted from the MINUIT documentation file.

//*-*

//*-*

//*-* command name ASSUMED to be in upper case

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



static MString cmd3 = "    ";



//*-*.......................................................................

//*-*

//*-*  Global HELP: Summary of all commands

//*-*  ====================================

//*-*

    if (comd[0] == '*') {

       fPrintf("   ==>List of MINUIT Interactive commands:");

       fPrintf(" CLEar     Reset all parameter names and values undefined");

       fPrintf(" CONtour   Make contour map of the user function");

       fPrintf(" EXIT      Exit from Interactive Minuit");

       fPrintf(" FIX       Cause parameter(s) to remain constant");

       fPrintf(" HESse     Calculate the Hessian or error matrix.");

       fPrintf(" IMPROVE   Search for a new minimum around current minimum");

       fPrintf(" MIGrad    Minimize by the method of Migrad");

       fPrintf(" MINImize  MIGRAD + SIMPLEX method if Migrad fails");

       fPrintf(" MINOs     Exact (non-linear) parameter error analysis");

       fPrintf(" MNContour Calculate one MINOS function contour");

       fPrintf(" PARameter Define or redefine new parameters and values");

       fPrintf(" RELease   Make previously FIXed parameters variable again");

       fPrintf(" REStore   Release last parameter fixed");

       fPrintf(" SAVe      Save current parameter values on a file");

       fPrintf(" SCAn      Scan the user function by varying parameters");

       fPrintf(" SEEk      Minimize by the method of Monte Carlo");

       fPrintf(" SET       Set various MINUIT constants or conditions");

       fPrintf(" SHOw      Show values of current constants or conditions");

       fPrintf(" SIMplex   Minimize by the method of Simplex");

	goto L99;

    }



    cmd3 = comd;

//*-* __________________________________________________________________

//*-*

//*-* --  Command CLEAR

//*-* --  =============

//*-*

    if (cmd3 == "CLE") {

       fPrintf(" ***>CLEAR");

       fPrintf(" Resets all parameter names and values to undefined.");

       fPrintf(" Must normally be followed by a PARameters command or ");

       fPrintf(" equivalent, in order to define parameter values.");

	goto L99;

    }

//*-* __________________________________________________________________

//*-* --

//*-* --  Command CONTOUR

//*-* --  ===============

//*-* .

    if (cmd3 == "CON") {

       fPrintf(" ***>CONTOUR <par1>  <par2>  [devs]  [ngrid]");

       fPrintf(" Instructs Minuit to trace contour lines of the user function");

       fPrintf(" with respect to the two parameters whose external numbers");

       fPrintf(" are <par1> and <par2>.");

       fPrintf(" Other variable parameters of the function, if any, will have");

       fPrintf(" their values fixed at the current values during the contour");

       fPrintf(" tracing. The optional parameter [devs] (default value 2.)");

       fPrintf(" gives the number of standard deviations in each parameter");

       fPrintf(" which should lie entirely within the plotting area.");

       fPrintf(" Optional parameter [ngrid] (default value 25 unless page");

       fPrintf(" size is too small) determines the resolution of the plot,");

       fPrintf(" i.e. the number of rows and columns of the grid at which the");

       fPrintf(" function will be evaluated. [See also MNContour.]");

	goto L99;

    }

//*-* __________________________________________________________________

//*-* --

//*-* --  Command END

//*-* --  ===========

//*-* .

    if (cmd3 == "END") {

       fPrintf(" ***>END");

       fPrintf(" Signals the end of a data block (i.e., the end of a fit),");

       fPrintf(" and implies that execution should continue, because another");

       fPrintf(" Data Block follows. A Data Block is a set of Minuit data");

       fPrintf(" consisting of");

       fPrintf("     (1) A Title,");

       fPrintf("     (2) One or more Parameter Definitions,");

       fPrintf("     (3) A blank line, and");

       fPrintf("     (4) A set of Minuit Commands.");

       fPrintf(" The END command is used when more than one Data Block is to");

       fPrintf(" be used with the same FCN function. It first causes Minuit");

       fPrintf(" to issue a CALL FCN with IFLAG=3, in order to allow FCN to");

       fPrintf(" perform any calculations associated with the final fitted");

       fPrintf(" parameter values, unless a CALL FCN 3 command has already");

       fPrintf(" been executed at the current FCN value.");

	goto L99;

    }

//*-* __________________________________________________________________

//*-* .

//*-* --

//*-* --  Command EXIT

//*-* --  ============

    if (cmd3 == "EXI") {

       fPrintf(" ***>EXIT");

       fPrintf(" Signals the end of execution.");

       fPrintf(" The EXIT command first causes Minuit to issue a CALL FCN");

       fPrintf(" with IFLAG=3, to allow FCN to perform any calculations");

       fPrintf(" associated with the final fitted parameter values, unless a");

       fPrintf(" CALL FCN 3 command has already been executed.");

	goto L99;

    }

//*-* __________________________________________________________________

//*-* --

//*-* --  Command FIX

//*-* --  ===========

//*-* .

    if (cmd3 == "FIX") {

       fPrintf(" ***>FIX} <parno> [parno] ... [parno]");

       fPrintf(" Causes parameter(s) <parno> to be removed from the list of");

       fPrintf(" variable parameters, and their value(s) will remain constant");

       fPrintf(" during subsequent minimizations, etc., until another command");

       fPrintf(" changes their value(s) or status.");

	goto L99;

    }

//*-* __________________________________________________________________

//*-* --

//*-* --  Command HESSE

//*-* --  =============

//*-* .

    if (cmd3 == "HES") {

       fPrintf(" ***>HESse  [maxcalls]");

       fPrintf(" Calculate, by finite differences, the Hessian or error matrix.");

       fPrintf("  That is, it calculates the full matrix of second derivatives");

       fPrintf(" of the function with respect to the currently variable");

       fPrintf(" parameters, and inverts it, printing out the resulting error");

       fPrintf(" matrix. The optional argument [maxcalls] specifies the");

       fPrintf(" (approximate) maximum number of function calls after which");

       fPrintf(" the calculation will be stopped.");

	goto L99;

    }

//*-* __________________________________________________________________

//*-* --

//*-* --  Command IMPROVE

//*-* --  ===============

//*-* .

    if (cmd3 == "IMP") {

       fPrintf(" ***>IMPROVE  [maxcalls]");

       fPrintf(" If a previous minimization has converged, and the current");

       fPrintf(" values of the parameters therefore correspond to a local");

       fPrintf(" minimum of the function, this command requests a search for");

       fPrintf(" additional distinct local minima.");

       fPrintf(" The optional argument [maxcalls] specifies the (approximate");

       fPrintf(" maximum number of function calls after which the calculation");

       fPrintf(" will be stopped.");

	goto L99;

    }

//*-* __________________________________________________________________

//*-* --

//*-* --  Command MIGRAD

//*-* --  ==============

//*-* .

    if (cmd3 == "MIG") {

       fPrintf(" ***>MIGrad  [maxcalls]  [tolerance]");

       fPrintf(" Causes minimization of the function by the method of Migrad,");

       fPrintf(" the most efficient and complete single method, recommended");

       fPrintf(" for general functions (see also MINImize).");

       fPrintf(" The minimization produces as a by-product the error matrix");

       fPrintf(" of the parameters, which is usually reliable unless warning");

       fPrintf(" messages are produced.");

       fPrintf(" The optional argument [maxcalls] specifies the (approximate)");

       fPrintf(" maximum number of function calls after which the calculation");

       fPrintf(" will be stopped even if it has not yet converged.");

       fPrintf(" The optional argument [tolerance] specifies required tolerance");

       fPrintf(" on the function value at the minimum.");

       fPrintf(" The default tolerance is 0.1, and the minimization will stop");

       fPrintf(" when the estimated vertical distance to the minimum (EDM) is");

       fPrintf(" less than 0.001*[tolerance]*UP (see [SET ERRordef]).");

	goto L99;

    }

//*-* __________________________________________________________________

//*-* --

//*-* --  Command MINIMIZE

//*-* --  ================

//*-* .

    if (comd == "MINI") {

       fPrintf(" ***>MINImize  [maxcalls] [tolerance]");

       fPrintf(" Causes minimization of the function by the method of Migrad,");

       fPrintf(" as does the MIGrad command, but switches to the SIMplex method");

       fPrintf(" if Migrad fails to converge. Arguments are as for MIGrad.");

       fPrintf(" Note that command requires four characters to be unambiguous.");

	goto L99;

    }

//*-* __________________________________________________________________

//*-* --

//*-* --  Command MINOS

//*-* --  =============

//*-* .

    if (comd == "MINO") {

       fPrintf(" ***>MINOs  [maxcalls]  [parno] [parno] ...");

       fPrintf(" Causes a Minos error analysis to be performed on the parameters");

       fPrintf(" whose numbers [parno] are specified. If none are specified,");

       fPrintf(" Minos errors are calculated for all variable parameters.");

       fPrintf(" Minos errors may be expensive to calculate, but are very");

       fPrintf(" reliable since they take account of non-linearities in the");

       fPrintf(" problem as well as parameter correlations, and are in general");

       fPrintf(" asymmetric.");

       fPrintf(" The optional argument [maxcalls] specifies the (approximate)");

       fPrintf(" maximum number of function calls per parameter requested,");

       fPrintf(" after which the calculation will stop for that parameter.");

	goto L99;

    }

//*-* __________________________________________________________________

//*-* --

//*-* --  Command MNCONTOUR

//*-* --  =================

//*-* .

    if (cmd3 == "MNC") {

       fPrintf(" ***>MNContour  <par1> <par2> [npts]");

       fPrintf(" Calculates one function contour of FCN with respect to");

       fPrintf(" parameters par1 and par2, with FCN minimized always with");

       fPrintf(" respect to all other NPAR-2 variable parameters (if any).");

       fPrintf(" Minuit will try to find npts points on the contour (default 20)");

       fPrintf(" If only two parameters are variable at the time, it is not");

       fPrintf(" necessary to specify their numbers. To calculate more than");

       fPrintf(" one contour, it is necessary to SET ERRordef to the appropriate");

       fPrintf(" value and issue the MNContour command for each contour.");

	goto L99;

    }

//*-* __________________________________________________________________

//*-* --

//*-* --  Command PARAMETER

//*-* --  =================

//*-* .

    if (cmd3 == "PAR") {

       fPrintf(" ***>PARameters");

       fPrintf(" followed by one or more parameter definitions.");

       fPrintf(" Parameter definitions are of the form:");

       fPrintf("   <number>  ''name''  <value>  <step>  [lolim] [uplim] ");

       fPrintf(" for example:");

       fPrintf("  3  ''K width''  1.2   0.1");

       fPrintf(" the last definition is followed by a blank line or a zero.");

	goto L99;

    }

//*-* __________________________________________________________________

//*-* --

//*-* --  Command RELEASE

//*-* --  ===============

//*-* .

    if (cmd3 == "REL") {

       fPrintf(" ***>RELease  <parno> [parno] ... [parno]");

       fPrintf(" If <parno> is the number of a previously variable parameter");

       fPrintf(" which has been fixed by a command: FIX <parno>, then that");

       fPrintf(" parameter will return to variable status.  Otherwise a warning");

       fPrintf(" message is printed and the command is ignored.");

       fPrintf(" Note that this command operates only on parameters which were");

       fPrintf(" at one time variable and have been FIXed. It cannot make");

       fPrintf(" constant parameters variable; that must be done by redefining");

       fPrintf(" the parameter with a PARameters command.");

	goto L99;

    }

//*-* __________________________________________________________________

//*-* --

//*-* --  Command RESTORE

//*-* --  ===============

//*-* .

    if (cmd3 == "RES") {

       fPrintf(" ***>REStore  [code]");

       fPrintf(" If no [code] is specified, this command restores all previously");

       fPrintf(" FIXed parameters to variable status. If [code]=1, then only");

       fPrintf(" the last parameter FIXed is restored to variable status.");

       fPrintf(" If code is neither zero nor one, the command is ignored.");

	goto L99;

    }

//*-* __________________________________________________________________

//*-* --

//*-* --  Command RETURN

//*-* --  ==============

//*-* .

    if (cmd3 == "RET") {

       fPrintf(" ***>RETURN");

       fPrintf(" Signals the end of a data block, and instructs Minuit to return");

       fPrintf(" to the program which called it. The RETurn command first");

       fPrintf(" causes Minuit to CALL FCN with IFLAG=3, in order to allow FCN");

       fPrintf(" to perform any calculations associated with the final fitted");

       fPrintf(" parameter values, unless a CALL FCN 3 command has already been");

       fPrintf(" executed at the current FCN value.");

	goto L99;

    }

//*-* __________________________________________________________________

//*-* --

//*-* --  Command SAVE

//*-* --  ============

//*-* .

    if (cmd3 == "SAV") {

       fPrintf(" ***>SAVe");

       fPrintf(" Causes the current parameter values to be saved on a file in");

       fPrintf(" such a format that they can be read in again as Minuit");

       fPrintf(" parameter definitions. If the covariance matrix exists, it is");

       fPrintf(" also output in such a format. The unit number is by default 7,");

       fPrintf(" or that specified by the user in his call to MINTIO or");

       fPrintf(" MNINIT. The user is responsible for opening the file previous");

       fPrintf(" to issuing the [SAVe] command (except where this can be done");

       fPrintf(" interactively).");

	goto L99;

    }

//*-* __________________________________________________________________

//*-* --

//*-* --  Command SCAN

//*-* --  ============

//*-* .

    if (cmd3 == "SCA") {

       fPrintf(" ***>SCAn  [parno]  [numpts] [from]  [to]");

       fPrintf(" Scans the value of the user function by varying parameter");

       fPrintf(" number [parno], leaving all other parameters fixed at the");

       fPrintf(" current value. If [parno] is not specified, all variable");

       fPrintf(" parameters are scanned in sequence.");

       fPrintf(" The number of points [numpts] in the scan is 40 by default,");

       fPrintf(" and cannot exceed 100. The range of the scan is by default");

       fPrintf(" 2 standard deviations on each side of the current best value,");

       fPrintf(" but can be specified as from [from] to [to].");

       fPrintf(" After each scan, if a new minimum is found, the best parameter");

       fPrintf(" values are retained as start values for future scans or");

       fPrintf(" minimizations. The curve resulting from each scan is plotted");

       fPrintf(" on the output unit in order to show the approximate behaviour");

       fPrintf(" of the function.");

       fPrintf(" This command is not intended for minimization, but is sometimes");

       fPrintf(" useful for debugging the user function or finding a");

       fPrintf(" reasonable starting point.");

	goto L99;

    }

//*-* __________________________________________________________________

//*-* --

//*-* --  Command SEEK

//*-* --  ============

//*-* .

    if (cmd3 == "SEE") {

       fPrintf(" ***>SEEk  [maxcalls]  [devs]");

       fPrintf(" Causes a Monte Carlo minimization of the function, by choosing");

       fPrintf(" random values of the variable parameters, chosen uniformly");

       fPrintf(" over a hypercube centered at the current best value.");

       fPrintf(" The region size is by default 3 standard deviations on each");

       fPrintf(" side, but can be changed by specifying the value of [devs].");

	goto L99;

    }

//*-* __________________________________________________________________

//*-* --

//*-* --  Command SET

//*-* --  ===========

//*-* .

    if (cmd3 == "SET") {

       fPrintf(" ***>SET <option_name>");

       fPrintf("  SET BATch");

       fPrintf("    Informs Minuit that it is running in batch mode.");



       fPrintf(" ");

       fPrintf("  SET EPSmachine  <accuracy>");

       fPrintf("    Informs Minuit that the relative floating point arithmetic");

       fPrintf("    precision is <accuracy>. Minuit determines the nominal");

       fPrintf("    precision itself, but the SET EPSmachine command can be");

       fPrintf("    used to override Minuit own determination, when the user");

       fPrintf("    knows that the FCN function value is not calculated to");

       fPrintf("    the nominal machine accuracy. Typical values of <accuracy>");

       fPrintf("    are between 10**-5 and 10**-14.");



       fPrintf(" ");

       fPrintf("  SET ERRordef  <up>");

       fPrintf("    Sets the value of UP (default value= 1.), defining parameter");

       fPrintf("    errors. Minuit defines parameter errors as the change");

       fPrintf("    in parameter value required to change the function value");

       fPrintf("    by UP. Normally, for chisquared fits UP=1, and for negative");

       fPrintf("    log likelihood, UP=0.5.");



       fPrintf(" ");

       fPrintf("   SET GRAdient  [force]");

       fPrintf("    Informs Minuit that the user function is prepared to");

       fPrintf("    calculate its own first derivatives and return their values");

       fPrintf("    in the array GRAD when IFLAG=2 (see specs of FCN).");

       fPrintf("    If [force] is not specified, Minuit will calculate");

       fPrintf("    the FCN derivatives by finite differences at the current");

       fPrintf("    point and compare with the user calculation at that point,");

       fPrintf("    accepting the user values only if they agree.");

       fPrintf("    If [force]=1, Minuit does not do its own derivative");

       fPrintf("    calculation, and uses the derivatives calculated in FCN.");



       fPrintf(" ");

       fPrintf("   SET INPut  [unitno]  [filename]");

       fPrintf("    Causes Minuit, in data-driven mode only, to read subsequent");

       fPrintf("    commands (or parameter definitions) from a different input");

       fPrintf("    file. If no [unitno] is specified, reading reverts to the");

       fPrintf("    previous input file, assuming that there was one.");

       fPrintf("    If [unitno] is specified, and that unit has not been opened,");

       fPrintf("    then Minuit attempts to open the file [filename]} if a");

       fPrintf("    name is specified. If running in interactive mode and");

       fPrintf("    [filename] is not specified and [unitno] is not opened,");

       fPrintf("    Minuit prompts the user to enter a file name.");

       fPrintf("    If the word REWIND is added to the command (note:no blanks");

       fPrintf("    between INPUT and REWIND), the file is rewound before");

       fPrintf("    reading. Note that this command is implemented in standard");

       fPrintf("    Fortran 77 and the results may depend on the  system;");

       fPrintf("    for example, if a filename is given under VM/CMS, it must");

       fPrintf("    be preceeded by a slash.");



       fPrintf(" ");

       fPrintf("   SET INTeractive");

       fPrintf("    Informs Minuit that it is running interactively.");



       fPrintf(" ");

       fPrintf("   SET LIMits  [parno]  [lolim]  [uplim]");

       fPrintf("    Allows the user to change the limits on one or all");

       fPrintf("    parameters. If no arguments are specified, all limits are");

       fPrintf("    removed from all parameters. If [parno] alone is specified,");

       fPrintf("    limits are removed from parameter [parno].");

       fPrintf("    If all arguments are specified, then parameter [parno] will");

       fPrintf("    be bounded between [lolim] and [uplim].");

       fPrintf("    Limits can be specified in either order, Minuit will take");

       fPrintf("    the smaller as [lolim] and the larger as [uplim].");

       fPrintf("    However, if [lolim] is equal to [uplim], an error condition");

       fPrintf("    results.");



       fPrintf(" ");

       fPrintf("   SET LINesperpage");

       fPrintf("     Sets the number of lines for one page of output.");

       fPrintf("     Default value is 24 for interactive mode");



       fPrintf(" ");

       fPrintf("   SET NOGradient");

       fPrintf("    The inverse of SET GRAdient, instructs Minuit not to");

       fPrintf("    use the first derivatives calculated by the user in FCN.");



       fPrintf(" ");

       fPrintf("   SET NOWarnings");

       fPrintf("    Supresses Minuit warning messages.");



       fPrintf(" ");

       fPrintf("   SET OUTputfile  <unitno>");

       fPrintf("    Instructs Minuit to write further output to unit <unitno>.");



       fPrintf(" ");

       fPrintf("   SET PAGethrow  <integer>");

       fPrintf("    Sets the carriage control character for ``new page'' to");

       fPrintf("    <integer>. Thus the value 1 produces a new page, and 0");

       fPrintf("    produces a blank line, on some devices (see TOPofpage)");





       fPrintf(" ");

       fPrintf("   SET PARameter  <parno>  <value>");

       fPrintf("    Sets the value of parameter <parno> to <value>.");

       fPrintf("    The parameter in question may be variable, fixed, or");

       fPrintf("    constant, but must be defined.");



       fPrintf(" ");

       fPrintf("   SET PRIntout  <level>");

       fPrintf("    Sets the print level, determining how much output will be");

       fPrintf("    produced. Allowed values and their meanings are displayed");

       fPrintf("    after a SHOw PRInt command, and are currently <level>=:");

       fPrintf("      [-1]  no output except from SHOW commands");

       fPrintf("       [0]  minimum output");

       fPrintf("       [1]  default value, normal output");

       fPrintf("       [2]  additional output giving intermediate results.");

       fPrintf("       [3]  maximum output, showing progress of minimizations.");

       fPrintf("    Note: See also the SET WARnings command.");



       fPrintf(" ");

       fPrintf("   SET RANdomgenerator  <seed>");

       fPrintf("    Sets the seed of the random number generator used in SEEk.");

       fPrintf("    This can be any integer between 10000 and 900000000, for");

       fPrintf("    example one which was output from a SHOw RANdom command of");

       fPrintf("    a previous run.");



       fPrintf(" ");

       fPrintf("   SET STRategy  <level>");

       fPrintf("    Sets the strategy to be used in calculating first and second");

       fPrintf("    derivatives and in certain minimization methods.");

       fPrintf("    In general, low values of <level> mean fewer function calls");

       fPrintf("    and high values mean more reliable minimization.");

       fPrintf("    Currently allowed values are 0, 1 (default), and 2.");



       fPrintf(" ");

       fPrintf("   SET TITle");

       fPrintf("    Informs Minuit that the next input line is to be considered");

       fPrintf("    the (new) title for this task or sub-task.  This is for");

       fPrintf("    the convenience of the user in reading his output.");



       fPrintf(" ");

       fPrintf("   SET WARnings");

       fPrintf("    Instructs Minuit to output warning messages when suspicious");

       fPrintf("    conditions arise which may indicate unreliable results.");

       fPrintf("    This is the default.");



       fPrintf(" ");

       fPrintf("    SET WIDthpage");

       fPrintf("    Informs Minuit of the output page width.");

       fPrintf("    Default values are 80 for interactive jobs");

	goto L99;

    }

//*-* __________________________________________________________________

//*-* --

//*-* --  Command SHOW

//*-* --  ============

//*-* .

    if (cmd3 == "SHO") {

       fPrintf(" ***>SHOw  <option_name>");

       fPrintf("  All SET XXXX commands have a corresponding SHOw XXXX command.");

       fPrintf("  In addition, the SHOw commands listed starting here have no");

       fPrintf("  corresponding SET command for obvious reasons.");



       fPrintf(" ");

       fPrintf("   SHOw CORrelations");

       fPrintf("    Calculates and prints the parameter correlations from the");

       fPrintf("    error matrix.");



       fPrintf(" ");

       fPrintf("   SHOw COVariance");

       fPrintf("    Prints the (external) covariance (error) matrix.");



       fPrintf(" ");

       fPrintf("   SHOw EIGenvalues");

       fPrintf("    Calculates and prints the eigenvalues of the covariance");

       fPrintf("    matrix.");



       fPrintf(" ");

       fPrintf("   SHOw FCNvalue");

       fPrintf("    Prints the current value of FCN.");

	goto L99;

    }

//*-* __________________________________________________________________

//*-* --

//*-* --  Command SIMPLEX

//*-* --  ===============

//*-* .

    if (cmd3 == "SIM") {

       fPrintf(" ***>SIMplex  [maxcalls]  [tolerance]");

       fPrintf(" Performs a function minimization using the simplex method of");

       fPrintf(" Nelder and Mead. Minimization terminates either when the");

       fPrintf(" function has been called (approximately) [maxcalls] times,");

       fPrintf(" or when the estimated vertical distance to minimum (EDM) is");

       fPrintf(" less than [tolerance].");

       fPrintf(" The default value of [tolerance] is 0.1*UP(see SET ERRordef).");

	goto L99;

    }

//*-* __________________________________________________________________

//*-* --

//*-* --  Command STANDARD

//*-* --  ================

//*-* .

    if (cmd3 == "STA") {

       fPrintf(" ***>STAndard");

	goto L99;

    }

//*-* __________________________________________________________________

//*-* --

//*-* --  Command STOP

//*-* --  ============

//*-* .

    if (cmd3 == "STO") {

       fPrintf(" ***>STOP");

       fPrintf(" Same as EXIT.");

	goto L99;

    }

//*-* __________________________________________________________________

//*-* --

//*-* --  Command TOPOFPAGE

//*-* --  =================

//*-* .

    if (cmd3 == "TOP") {

       fPrintf(" ***>TOPofpage");

       fPrintf(" Causes Minuit to write the character specified in a");

       fPrintf(" SET PAGethrow command (default = 1) to column 1 of the output");

       fPrintf(" file, which may or may not position your output medium to");

       fPrintf(" the top of a page depending on the device and system.");

	goto L99;

    }

//*-* __________________________________________________________________

       fPrintf(" Unknown MINUIT command. Type HELP for list of commands.");



L99:

    return;

} /* mnhelp_ */



//______________________________________________________________________________

void Midnight::mnhess()

{

//*-*-*-*-*-*Calculates the full second-derivative matrix of FCN*-*-*-*-*-*-*-*

//*-*        ===================================================

//*-*        by taking finite differences. When calculating diagonal

//*-*        elements, it may iterate so that step size is nearly that

//*-*        which gives function change= UP/10. The first derivatives

//*-*        of course come as a free side effect, but with a smaller

//*-*        step size in order to obtain a known accuracy.

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    /* Local variables */

    static MDouble dmin_, dxdi, elem, wint, tlrg2, d, dlast, ztemp, g2bfor;

    static MDouble  yy[kMAXDIM], df, aimsag, fs1, tlrstp, fs2, stpinm, g2i, sag, xtf, xti, xtj;

    static MInt icyc, ncyc, ndex, idrv, iext, npar2, i, j, ifail, npard, nparx, id, multpy;

    static MBool ldebug;



    ldebug = fIdbg[3] >= 1;

    if (fAmin == fUndefi) {

	mnamin();

    }

    if (fIstrat <= 0) {

	ncyc   = 3;

	tlrstp = .5;

	tlrg2  = .1;

    } else if (fIstrat == 1) {

	ncyc   = 5;

	tlrstp = .3;

	tlrg2  = .05;

    } else {

	ncyc   = 7;

	tlrstp = .1;

	tlrg2  = .02;

    }

    if (fISW[4] >= 2 || ldebug) {

	fPrintf("   START COVARIANCE MATRIX CALCULATION.");

    }

    fCfrom  = "HESSE   ";

    fNfcnfr = fNfcn;

    fCstatu = "OK        ";

    npard   = fNpar;

//*-*-                make sure starting at the right place

    mninex(fX);

    nparx = fNpar;

    (*fFCN)(nparx, fGin, fs1, fU, 4);    ++fNfcn;

    if (fs1 != fAmin) {

        df    = fAmin - fs1;

	mnwarn("D", "MNHESS", Form("function value differs from AMIN by %g",df));

    }

    fAmin = fs1;

    if (ldebug) {

	fPrintf(" PAR D   GSTEP           D          G2         GRD         SAG    ");

    }

//*-*-                                       . . . . . . diagonal elements .



//*-*-        ISW(2) = 1 if approx, 2 if not posdef, 3 if ok

//*-*-        AIMSAG is the sagitta we are aiming for in second deriv calc.



    aimsag = TMath_Sqrt(fEpsma2)*(TMath_Abs(fAmin) + fUp);

//*-*-        Zero the second derivative matrix

    npar2 = fNpar*(fNpar + 1) / 2;

    for (i = 1; i <= npar2; ++i) { fVhmat[i-1] = 0; }



//*-*-        Loop over variable parameters for second derivatives

    idrv = 2;

    for (id = 1; id <= npard; ++id) {

	i = id + fNpar - npard;

	iext = fNexofi[i-1];

	if (fG2[i-1] == 0) {

	    mnwarn("W", "HESSE", Form("Second derivative enters zero, param %d",iext));

	    wint = fWerr[i-1];

	    if (fNvarl[iext-1] > 1) {

		mndxdi(fX[i-1], i-1, dxdi);

		if (TMath_Abs(dxdi) < .001) wint = .01;

		else                          wint /= TMath_Abs(dxdi);

	    }

	    fG2[i-1] = fUp / (wint*wint);

	}

	xtf   = fX[i-1];

	dmin_ = fEpsma2*8*TMath_Abs(xtf);



//*-*-                              find step which gives sagitta = AIMSAG

	d = TMath_Abs(fGstep[i-1]);

        int skip50 = 0;

	for (icyc = 1; icyc <= ncyc; ++icyc) {

//*-*-                              loop here only if SAG=0

	    for (multpy = 1; multpy <= 5; ++multpy) {

//*-*-          take two steps

		fX[i-1] = xtf + d;

		mninex(fX);

		nparx = fNpar;

		(*fFCN)(nparx, fGin, fs1, fU, 4);    ++fNfcn;

		fX[i-1] = xtf - d;

		mninex(fX);

		(*fFCN)(nparx, fGin, fs2, fU, 4);    ++fNfcn;

		fX[i-1] = xtf;

		sag = (fs1 + fs2 - fAmin*2)*.5;

		if (sag != 0) goto L30;

		if (fGstep[i-1] < 0) {

		    if (d >= .5) goto L26;

		    d *= 10;

		    if (d > .5) 	d = .51;

		    continue;

		}

		d *= 10;

	    }

L26:

	    mnwarn("W", "HESSE", Form("Second derivative zero for parameter%d",iext));

	    goto L390;

//*-*-                            SAG is not zero

L30:

	    g2bfor    = fG2[i-1];

	    fG2[i-1]  = sag*2 / (d*d);

	    fGrd[i-1] = (fs1 - fs2) / (d*2);

	    if (ldebug) {

		fPrintf("%4d%2d%12.5g%12.5g%12.5g%12.5g%12.5g%12.5g",i,idrv,fGstep[i-1],fG2[i-1],fGrd[i-1],sag);

	    }

	    if (fGstep[i-1] > 0) fGstep[i-1] =  TMath_Abs(d);

	    else                 fGstep[i-1] = -TMath_Abs(d);

	    fDirin[i-1] = d;

	    yy[i-1]     = fs1;

	    dlast       = d;

	    d           = TMath_Sqrt(aimsag*2 / TMath_Abs(fG2[i-1]));

//*-*-        if parameter has limits, max int step size = 0.5

	    stpinm = .5;

	    if (fGstep[i-1] < 0) d = TMath_Min(d,stpinm);

	    if (d < dmin_) d = dmin_;

//*-*-          see if converged

	    if (TMath_Abs((d - dlast) / d) < tlrstp ||

	        TMath_Abs((fG2[i-1] - g2bfor) / fG2[i-1]) < tlrg2) {

               skip50 = 1;

               break;

            }

	    d = TMath_Min(d,dlast*102);

	    d = TMath_Max(d,dlast*.1);

	}

//*-*-                      end of step size loop

        if (!skip50)

	   mnwarn("D", "MNHESS", Form("Second Deriv. SAG,AIM= %d%g%g",iext,sag,aimsag));



	ndex = i*(i + 1) / 2;

	fVhmat[ndex-1] = fG2[i-1];

    }

//*-*-                             end of diagonal second derivative loop

    mninex(fX);

//*-*-                                    refine the first derivatives

    if (fIstrat > 0) mnhes1();

    fISW[1] = 3;

    fDcovar = 0;

//*-*-                                       . . . .  off-diagonal elements



    if (fNpar == 1) goto L214;

    for (i = 1; i <= fNpar; ++i) {

	for (j = 1; j <= i-1; ++j) {

	    xti     = fX[i-1];

	    xtj     = fX[j-1];

	    fX[i-1] = xti + fDirin[i-1];

	    fX[j-1] = xtj + fDirin[j-1];

	    mninex(fX);

	    (*fFCN)(nparx, fGin, fs1, fU, 4);	    ++fNfcn;

	    fX[i-1] = xti;

	    fX[j-1] = xtj;

	    elem = (fs1 + fAmin - yy[i-1] - yy[j-1]) / (

		    fDirin[i-1]*fDirin[j-1]);

	    ndex = i*(i-1) / 2 + j;

	    fVhmat[ndex-1] = elem;

	}

    }

L214:

    mninex(fX);

//*-*-                 verify matrix positive-definite

    mnpsdf();

    for (i = 1; i <= fNpar; ++i) {

	for (j = 1; j <= i; ++j) {

	    ndex = i*(i-1) / 2 + j;

	    fP[i + j*fMaxpar - fMaxpar-1] = fVhmat[ndex-1];

	    fP[j + i*fMaxpar - fMaxpar-1] = fP[i + j*fMaxpar - fMaxpar-1];

	}

    }

    mnvert(fP, fMaxint, fMaxint, fNpar, ifail);

    if (ifail > 0) {

	mnwarn("W", "HESSE", "Matrix inversion fails.");

	goto L390;

    }

//*-*-                                       . . . . . . .  calculate  e d m

    fEDM = 0;



    for (i = 1; i <= fNpar; ++i) {

//*-*-                             off-diagonal elements

	ndex = i*(i-1) / 2;

	for (j = 1; j <= i-1; ++j) {

	    ++ndex;

	    ztemp = fP[i + j*fMaxpar - fMaxpar-1]*2;

	    fEDM += fGrd[i-1]*ztemp*fGrd[j-1];

	    fVhmat[ndex-1] = ztemp;

	}

//*-*-                             diagonal elements

	++ndex;

	fVhmat[ndex-1] = fP[i + i*fMaxpar - fMaxpar-1]*2;

	fEDM += fP[i + i*fMaxpar - fMaxpar-1]*(fGrd[i-1]*fGrd[i-1]);

    }

    if (fISW[4] >= 1 && fISW[1] == 3 && fItaur == 0) {

	fPrintf(" COVARIANCE MATRIX CALCULATED SUCCESSFULLY");

    }

    goto L900;

//*-*-                             failure to invert 2nd deriv matrix

L390:

    fISW[1] = 1;

    fDcovar = 1;

    fCstatu = "FAILED    ";

    if (fISW[4] >= 0) {

	fPrintf("  MNHESS FAILS AND WILL RETURN DIAGONAL MATRIX. ");

    }

    for (i = 1; i <= fNpar; ++i) {

	ndex = i*(i-1) / 2;

	for (j = 1; j <= i-1; ++j) {

	    ++ndex;

	    fVhmat[ndex-1] = 0;

	}

	++ndex;

	g2i = fG2[i-1];

	if (g2i <= 0) g2i = 1;

	fVhmat[ndex-1] = 2 / g2i;

    }

L900:

    return;

} /* mnhess_ */



//______________________________________________________________________________

void Midnight::mnhes1()

{

//*-*-*-*Calculate first derivatives (GRD) and uncertainties (DGRD)*-*-*-*-*-*

//*-*    ==========================================================

//*-*         and appropriate step sizes GSTEP

//*-*      Called from MNHESS and MNGRAD

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    /* Local variables */

    static MDouble dmin_, d, dfmin, dgmin, change, chgold, grdold, epspri;

    static MDouble fs1, optstp, fs2, grdnew, sag, xtf;

    static MInt icyc, ncyc, idrv, i, nparx;

    static MBool ldebug;



    ldebug = fIdbg[5] >= 1;

    if (fIstrat <= 0) ncyc = 1;

    if (fIstrat == 1) ncyc = 2;

    if (fIstrat > 1)  ncyc = 6;

    idrv = 1;

    nparx = fNpar;

    dfmin = fEpsma2*4*(TMath_Abs(fAmin) + fUp);

//*-*-                                    main loop over parameters

    for (i = 1; i <= fNpar; ++i) {

	xtf    = fX[i-1];

	dmin_  = fEpsma2*4*TMath_Abs(xtf);

	epspri = fEpsma2 + TMath_Abs(fGrd[i-1]*fEpsma2);

	optstp = TMath_Sqrt(dfmin / (TMath_Abs(fG2[i-1]) + epspri));

	d = TMath_Abs(fGstep[i-1])*.2;

	if (d > optstp) d = optstp;

	if (d < dmin_)  d = dmin_;

	chgold = 1e4;

//*-*-                                      iterate reducing step size

	for (icyc = 1; icyc <= ncyc; ++icyc) {

	    fX[i-1] = xtf + d;

	    mninex(fX);

	    (*fFCN)(nparx, fGin, fs1, fU, 4);	    ++fNfcn;

	    fX[i-1] = xtf - d;

	    mninex(fX);

	    (*fFCN)(nparx, fGin, fs2, fU, 4);	    ++fNfcn;

	    fX[i-1] = xtf;

//*-*-                                      check if step sizes appropriate

	    sag    = (fs1 + fs2 - fAmin*2)*.5;

	    grdold = fGrd[i-1];

	    grdnew = (fs1 - fs2) / (d*2);

	    dgmin  = fEpsmac*(TMath_Abs(fs1) + TMath_Abs(fs2)) / d;

	    if (ldebug) {

		fPrintf("%4d%2d%12.5g%12.5g%12.5g%12.5g%12.5g%12.5g",i,idrv,fGstep[i-1],d,fG2[i-1],grdnew,sag);

	    }

	    if (grdnew == 0) goto L60;

	    change = TMath_Abs((grdold - grdnew) / grdnew);

	    if (change > chgold && icyc > 1) goto L60;

	    chgold    = change;

	    fGrd[i-1] = grdnew;

	    if (fGstep[i-1] > 0) fGstep[i-1] =  TMath_Abs(d);

	    else                 fGstep[i-1] = -TMath_Abs(d);

//*-*-                 decrease step until first derivative changes by <5%

	    if (change < .05) goto L60;

	    if (TMath_Abs(grdold - grdnew) < dgmin) goto L60;

	    if (d < dmin_) {

		mnwarn("D", "MNHES1", "Step size too small for 1st drv.");

		goto L60;

	    }

	    d *= .2;

	}

//*-*-                                      loop satisfied = too many iter

	mnwarn("D", "MNHES1", Form("Too many iterations on D1.%g%g",grdold,grdnew));

L60:

	fDgrd[i-1] = TMath_Max(dgmin,TMath_Abs(grdold - grdnew));

    }

//*-*-                                       end of first deriv. loop

    mninex(fX);

} /* mnhes1_ */



//______________________________________________________________________________

void Midnight::mnimpr()

{

//*-*-*-*-*-*-*Attempts to improve on a good local minimum*-*-*-*-*-*-*-*-*-*

//*-*          ===========================================

//*-*        Attempts to improve on a good local minimum by finding a

//*-*        better one.   The quadratic part of FCN is removed by MNCALF

//*-*        and this transformed function is minimized using the simplex

//*-*        method from several random starting points.

//*-*        ref. -- Goldstein and Price, Math.Comp. 25, 569 (1971)

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    /* Initialized data */



    static MDouble rnum = 0;



    /* Local variables */

    static MDouble amax, dsav[kMAXDIM], y[kMAXDIM+1], ycalf, ystar, ystst;

    static MDouble pb, ep, wg, xi, sigsav, reg, sig2;

    static MInt npfn, ndex, loop, i, j, ifail, iseed;

    static MInt jhold, nloop, nparx, nparp1, jh, jl;



    if (fNpar <= 0) return;

    if (fAmin == fUndefi) mnamin();

    fCstatu = "UNCHANGED ";

    fItaur  = 1;

    fEpsi   = fUp*.1;

    npfn    = fNfcn;

    nloop   = MInt(fWord7[1]);

    if (nloop <= 0) nloop = fNpar + 4;

    nparx  = fNpar;

    nparp1 = fNpar + 1;

    wg = 1 / MDouble(fNpar);

    sigsav = fEDM;

    fApsi  = fAmin;

    for (i = 1; i <= fNpar; ++i) {

	fXt[i-1]  = fX[i-1];

	dsav[i-1] = fWerr[i-1];

	for (j = 1; j <= i; ++j) {

	    ndex = i*(i-1) / 2 + j;

	    fP[i + j*fMaxpar - fMaxpar-1] = fVhmat[ndex-1];

	    fP[j + i*fMaxpar - fMaxpar-1] = fP[i + j*fMaxpar - fMaxpar-1];

	}

    }

    mnvert(fP, fMaxint, fMaxint, fNpar, ifail);

    if (ifail >= 1) goto L280;

//*-*-              Save inverted matrix in VT

    for (i = 1; i <= fNpar; ++i) {

	ndex = i*(i-1) / 2;

	for (j = 1; j <= i; ++j) {

	    ++ndex;

	    fVthmat[ndex-1] = fP[i + j*fMaxpar - fMaxpar-1];

	}

    }

    loop = 0;



L20:

    for (i = 1; i <= fNpar; ++i) {

	fDirin[i-1] = dsav[i-1]*2;

	mnrn15(rnum, iseed);

	fX[i-1] = fXt[i-1] + fDirin[i-1]*2*(rnum - .5);

    }

    ++loop;

    reg = 2;

    if (fISW[4] >= 0) {

	fPrintf("START ATTEMPT NO.%2d TO FIND NEW MINIMUM",loop);

    }

L30:

    mncalf(fX, ycalf);

    fAmin = ycalf;

//*-*-                                       . . . . set up  random simplex

    jl = nparp1;

    jh = nparp1;

    y[nparp1-1] = fAmin;

    amax = fAmin;

    for (i = 1; i <= fNpar; ++i) {

	xi = fX[i-1];

	mnrn15(rnum, iseed);

	fX[i-1] = xi - fDirin[i-1]*(rnum - .5);

	mncalf(fX, ycalf);

	y[i-1] = ycalf;

	if (y[i-1] < fAmin) {

	    fAmin = y[i-1];

	    jl    = i;

	} else if (y[i-1] > amax) {

	    amax = y[i-1];

	    jh   = i;

	}

	for (j = 1; j <= fNpar; ++j) { fP[j + i*fMaxpar - fMaxpar-1] = fX[j-1]; }

	fP[i + nparp1*fMaxpar - fMaxpar-1] = xi;

	fX[i-1] = xi;

    }



    fEDM = fAmin;

    sig2 = fEDM;

//*-*-                                       . . . . . . .  start main loop

L50:

    if (fAmin < 0)   goto L95;

    if (fISW[1] <= 2) goto L280;

    ep = fAmin*.1;

    if (sig2 < ep && fEDM < ep) goto L100;

    sig2 = fEDM;

    if (fNfcn - npfn > fNfcnmx) goto L300;

//*-*-        calculate new point * by reflection

    for (i = 1; i <= fNpar; ++i) {

	pb = 0;

	for (j = 1; j <= nparp1; ++j) { pb += wg*fP[i + j*fMaxpar - fMaxpar-1]; }

	fPbar[i-1]  = pb - wg*fP[i + jh*fMaxpar - fMaxpar-1];

	fPstar[i-1] = fPbar[i-1]*2 - fP[i + jh*fMaxpar - fMaxpar-1]*1;

    }

    mncalf(fPstar, ycalf);

    ystar = ycalf;

    if (ystar >= fAmin) goto L70;

//*-*-        point * better than jl, calculate new point **

    for (i = 1; i <= fNpar; ++i) {

	fPstst[i-1] = fPstar[i-1]*2 + fPbar[i- 1]*-1;

    }

    mncalf(fPstst, ycalf);

    ystst = ycalf;

    if (ystst < y[jl-1]) goto L67;

    mnrazz(ystar, fPstar, y, jh, jl);

    goto L50;

L67:

    mnrazz(ystst, fPstst, y, jh, jl);

    goto L50;

//*-*-        point * is not as good as jl

L70:

    if (ystar >= y[jh-1]) goto L73;

    jhold = jh;

    mnrazz(ystar, fPstar, y, jh, jl);

    if (jhold != jh) goto L50;

//*-*-        calculate new point **

L73:

    for (i = 1; i <= fNpar; ++i) {

	fPstst[i-1] = fP[i + jh*fMaxpar - fMaxpar-1]*.5 + fPbar[i-1]*.5;

    }

    mncalf(fPstst, ycalf);

    ystst = ycalf;

    if (ystst > y[jh-1]) goto L30;

//*-*-    point ** is better than jh

    if (ystst < fAmin) goto L67;

    mnrazz(ystst, fPstst, y, jh, jl);

    goto L50;

//*-*-                                       . . . . . .  end main loop

L95:

    if (fISW[4] >= 0) {

	fPrintf(" AN IMPROVEMENT ON THE PREVIOUS MINIMUM HAS BEEN FOUND");

    }

    reg = .1;

//*-*-                                       . . . . . ask if point is new

L100:

    mninex(fX);

    (*fFCN)(nparx, fGin, fAmin, fU, 4);    ++fNfcn;

    for (i = 1; i <= fNpar; ++i) {

	fDirin[i-1] = reg*dsav[i-1];

	if (TMath_Abs(fX[i-1] - fXt[i-1]) > fDirin[i-1])     goto L150;

    }

    goto L230;

L150:

    fNfcnmx = fNfcnmx + npfn - fNfcn;

    npfn    = fNfcn;

    mnsimp();

    if (fAmin >= fApsi) goto L325;

    for (i = 1; i <= fNpar; ++i) {

	fDirin[i-1] = dsav[i-1]*.1;

	if (TMath_Abs(fX[i-1] - fXt[i-1]) > fDirin[i-1])     goto L250;

    }

L230:

    if (fAmin < fApsi) 	goto L350;

    goto L325;

/*                                        . . . . . . truly new minimum */

L250:

    fLnewmn = kTRUE;

    if (fISW[1] >= 1) {

	fISW[1] = 1;

	fDcovar = TMath_Max(fDcovar,.5);

    } else fDcovar = 1;

    fItaur  = 0;

    fNfcnmx = fNfcnmx + npfn - fNfcn;

    fCstatu = "NEW MINIMU";

    if (fISW[4] >= 0) {

	fPrintf(" IMPROVE HAS FOUND A TRULY NEW MINIMUM");

        fPrintf(" *************************************");

    }

    return;

//*-*-                                       . . . return to previous region

L280:

    if (fISW[4] > 0) {

	fPrintf(" COVARIANCE MATRIX WAS NOT POSITIVE-DEFINITE");

    }

    goto L325;

L300:

    fISW[0] = 1;

L325:

    for (i = 1; i <= fNpar; ++i) {

	fDirin[i-1] = dsav[i-1]*.01;

	fX[i-1]     = fXt[i-1];

    }

    fAmin = fApsi;

    fEDM  = sigsav;

L350:

    mninex(fX);

    if (fISW[4] > 0) {

	fPrintf(" IMPROVE HAS RETURNED TO REGION OF ORIGINAL MINIMUM");

    }

    fCstatu = "UNCHANGED ";

    mnrset(0);

    if (fISW[1] < 2) goto L380;

    if (loop < nloop && fISW[0] < 1) goto L20;

L380:

    mnprin(5, fAmin);

    fItaur = 0;

} /* mnimpr_ */



//______________________________________________________________________________

void Midnight::mninex(MDouble *pint)

{

//*-*-*-*-*Transforms from internal coordinates (PINT) to external (U)*-*-*-*

//*-*      ===========================================================

//*-*        The minimizing routines which work in

//*-*        internal coordinates call this routine before calling FCN.

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    MInt i, j;



    for (j = 1; j <= fNpar; ++j) {

	i = fNexofi[j-1];

	if (fNvarl[i-1] == 1) {

	    fU[i-1] = pint[j-1];

	} else {

	    fU[i-1] = fAlim[i-1] + (TMath_Sin(pint[j-1]) + 1)*.5*(fBlim[i-1] - fAlim[i-1]);

	}

    }

} /* mninex_ */



//______________________________________________________________________________

void Midnight::mninit(MInt i1, MInt i2, MInt i3)

{

//*-*-*-*-*-*Main initialization member function for MINUIT*-*-*-*-*-*-*-*-*

//*-*        ==============================================

//*-*     It initializes some constants

//*-*                (including the logical I/O unit nos.),

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    /* Local variables */

    static MDouble piby2, epsp1, epsbak, epstry, distnn;

    static MInt i, idb;



//*-*-           I/O unit numbers

    fIsysrd = i1;

    fIsyswr = i2;

    fIstkwr[0] = fIsyswr;

    fNstkwr = 1;

    fIsyssa = i3;

    fNstkrd = 0;

//*-*-              version identifier

    fCvrsn  = "95.03++ ";

//*-*-              some CONSTANT

    fMaxint = fMaxpar;

    fMaxext = 2*fMaxpar;

    fUndefi = -54321;

    fBigedm = 123456;

    fCundef = ")UNDEFINED";

    fCovmes[0] = "NO ERROR MATRIX       ";

    fCovmes[1] = "ERR MATRIX APPROXIMATE";

    fCovmes[2] = "ERR MATRIX NOT POS-DEF";

    fCovmes[3] = "ERROR MATRIX ACCURATE ";

//*-*-               some starting values

    fNblock = 0;

    fIcomnd = 0;

    fCtitl  = fCundef;

    fCfrom  = "INPUT   ";

    fNfcn   = 0;

    fNfcnfr = fNfcn;

    fCstatu = "INITIALIZE";

    fISW[2] = 0;

    fISW[3] = 0;

    fISW[4] = 1;

//*-*-        ISW(6)=0 for batch jobs,  =1 for interactive jobs

//*-*-                     =-1 for originally interactive temporarily batch



    fISW[5] = 0;

//    if (intrac(&dummy)) fISW[5] = 1;

//*-*-       DEBUG options set to default values

    for (idb = 0; idb <= 10; ++idb) { fIdbg[idb] = 0; }

    fLrepor = kFALSE;

    fLwarn  = kTRUE;

    fLimset = kFALSE;

    fLnewmn = kFALSE;

    fIstrat = 1;

    fItaur  = 0;

//*-*-       default page dimensions and 'new page' carriage control integer

    fNpagwd = 120;

    fNpagln = 56;

    fNewpag = 1;

    if (fISW[5] > 0) {

	fNpagwd = 80;

	fNpagln = 30;

	fNewpag = 0;

    }

    fUp = 1;

    fUpdflt = fUp;

//*-*-                  determine machine accuracy epsmac

    epstry = .5;

    for (i = 1; i <= 100; ++i) {

	epstry *= .5;

	epsp1 = epstry + 1;

	mntiny(epsp1, epsbak);

	if (epsbak < epstry) goto L35;

    }

    epstry = 1e-7;

    fEpsmac = epstry*4;

    fPrintf(" MNINIT UNABLE TO DETERMINE ARITHMETIC PRECISION. WILL ASSUME:%g",fEpsmac);

L35:

    fEpsmac = epstry*8;

    fEpsma2 = TMath_Sqrt(fEpsmac)*2;

//*-*-                the vlims are a non-negligible distance from pi/2

//*-*-        used by MNPINT to set variables "near" the physical limits

    piby2   = TMath_ATan(1.)*2;

    distnn  = TMath_Sqrt(fEpsma2)*8;

    fVlimhi =  piby2 - distnn;

    fVlimlo = -piby2 + distnn;

    mncler();

//    fPrintf("  MINUIT RELEASE %s INITIALIZED.   DIMENSIONS 100/50  EPSMAC=%g",(const char*)fCvrsn,fEpsmac);

} /* mninit_ */



//______________________________________________________________________________

void Midnight::mnlims()

{

//*-*-*-*Interprets the SET LIM command, to reset the parameter limits*-*-*-*

//*-*    =============================================================

//*-*       Called from MNSET

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    /* Local variables */

    static MDouble dxdi, snew;

    static MInt kint, i2, newcod, ifx, inu;



    fCfrom  = "SET LIM ";

    fNfcnfr = fNfcn;

    fCstatu = "NO CHANGE ";

    i2 = MInt(fWord7[0]);

    if (i2 > fMaxext || i2 < 0) goto L900;

    if (i2 > 0) goto L30;

//*-*-                                    set limits on all parameters

    newcod = 4;

    if (fWord7[1] == fWord7[2]) newcod = 1;

    for (inu = 1; inu <= fNu; ++inu) {

	if (fNvarl[inu-1] <= 0) continue;

	if (fNvarl[inu-1] == 1 && newcod == 1) continue;

	kint = fNiofex[inu-1];

//*-*-            see if parameter has been fixed

	if (kint <= 0) {

	    if (fISW[4] >= 0) {

		fPrintf("           LIMITS NOT CHANGED FOR FIXED PARAMETER:%4d",inu);

	    }

	    continue;

	}

	if (newcod == 1) {

//*-*-           remove limits from parameter

	    if (fISW[4] > 0) {

                fPrintf(" LIMITS REMOVED FROM PARAMETER  :%3d",inu);

	    }

	    fCstatu = "NEW LIMITS";

	    mndxdi(fX[kint-1], kint-1, dxdi);

	    snew           = fGstep[kint-1]*dxdi;

	    fGstep[kint-1] = TMath_Abs(snew);

	    fNvarl[inu-1]  = 1;

	} else {

//*-*-            put limits on parameter

	    fAlim[inu-1] = TMath_Min(fWord7[1],fWord7[2]);

	    fBlim[inu-1] = TMath_Max(fWord7[1],fWord7[2]);

	    if (fISW[4] > 0) {

		fPrintf(" PARAMETER %3d LIMITS SET TO  %15.5g%15.5g",inu,fAlim[inu-1],fBlim[inu-1]);

	    }

	    fNvarl[inu-1]  = 4;

	    fCstatu        = "NEW LIMITS";

	    fGstep[kint-1] = -.1;

	}

    }

    goto L900;

//*-*-                                      set limits on one parameter

L30:

    if (fNvarl[i2-1] <= 0) {

        fPrintf(" PARAMETER %3d IS NOT VARIABLE.", i2);

	goto L900;

    }

    kint = fNiofex[i2-1];

//*-*-                                      see if parameter was fixed

    if (kint == 0) {

        fPrintf(" REQUEST TO CHANGE LIMITS ON FIXED PARAMETER:%3d",i2);

	for (ifx = 1; ifx <= fNpfix; ++ifx) {

	    if (i2 == fIpfix[ifx-1]) goto L92;

	}

	fPrintf(" MINUIT BUG IN MNLIMS. SEE F. JAMES");

L92:

	;

    }

    if (fWord7[1] != fWord7[2]) goto L235;

//*-*-                                      remove limits

    if (fNvarl[i2-1] != 1) {

	if (fISW[4] > 0) {

            fPrintf(" LIMITS REMOVED FROM PARAMETER  %2d",i2);

	}

	fCstatu = "NEW LIMITS";

	if (kint <= 0) {

	    fGsteps[ifx-1] = TMath_Abs(fGsteps[ifx-1]);

	} else {

	    mndxdi(fX[kint-1], kint-1, dxdi);

	    if (TMath_Abs(dxdi) < .01) dxdi = .01;

	    fGstep[kint-1] = TMath_Abs(fGstep[kint-1]*dxdi);

	    fGrd[kint-1]  *= dxdi;

	}

	fNvarl[i2-1] = 1;

    } else {

	fPrintf(" NO LIMITS SPECIFIED.  PARAMETER %3d IS ALREADY UNLIMITED.  NO CHANGE.",i2);

    }

    goto L900;

//*-*-                                       put on limits

L235:

    fAlim[i2-1]  = TMath_Min(fWord7[1],fWord7[2]);

    fBlim[i2-1]  = TMath_Max(fWord7[1],fWord7[2]);

    fNvarl[i2-1] = 4;

    if (fISW[4] > 0) {

	fPrintf(" PARAMETER %3d LIMITS SET TO  %15.5g%15.5g",i2,fAlim[i2-1],fBlim[i2-1]);

    }

    fCstatu = "NEW LIMITS";

    if (kint <= 0) fGsteps[ifx-1] = -.1;

    else           fGstep[kint-1] = -.1;



L900:

    if (fCstatu != "NO CHANGE ") {

	mnexin(fX);

	mnrset(1);

    }

} /* mnlims_ */



//______________________________________________________________________________

void Midnight::mnline(MDouble *start, MDouble fstart, MDouble *step, MDouble slope, MDouble toler)

{

//*-*-*-*-*-*-*-*-*-*Perform a line search from position START*-*-*-*-*-*-*-*

//*-*                =========================================

//*-*        along direction STEP, where the length of vector STEP

//*-*                   gives the expected position of minimum.

//*-*        FSTART is value of function at START

//*-*        SLOPE (if non-zero) is df/dx along STEP at START

//*-*        TOLER is initial tolerance of minimum in direction STEP

//*-*

//*-* SLAMBG and ALPHA control the maximum individual steps allowed.

//*-* The first step is always =1. The max length of second step is SLAMBG.

//*-* The max size of subsequent steps is the maximum previous successful

//*-*   step multiplied by ALPHA + the size of most recent successful step,

//*-*   but cannot be smaller than SLAMBG.

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    const MString charal = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";



    /* Local variables */

    static MDouble xpq[12], ypq[12], slam, sdev, coeff[3], denom, flast;

    static MDouble fvals[3], xvals[3], f1, fvmin, xvmin, ratio, f2, f3, fvmax;

    static MDouble toler8, toler9, overal, undral, slamin, slamax, slopem;

    static MInt i, nparx, nvmax, nxypt, kk, ipt;

    static MBool ldebug;

    MString cmess, chpq[12];

    int     l65, l70, l80;



    /* Function Body */

    l65 = 0; l70 = 0; l80 = 0;

    ldebug = fIdbg[1] >= 1;

//*-*-                 starting values for overall limits on total step SLAM

    overal = 1e3;

    undral = -100;

//*-*-                             debug check if start is ok

    if (ldebug) {

	mninex(&start[0]);

	(*fFCN)(nparx, fGin, f1, fU, 4);	++fNfcn;

	if (f1 != fstart) {

	    fPrintf(" MNLINE start point not consistent, F values, parameters=");

	    for (kk = 1; kk <= fNpar; ++kk) {

                fPrintf("  %14.5e",fX[kk-1]);

	    }

	}

    }

//*-*-                                     . set up linear search along STEP

    fvmin   = fstart;

    xvmin   = 0;

    nxypt   = 1;

    chpq[0] = charal[0];

    xpq[0]  = 0;

    ypq[0]  = fstart;

//*-*-              SLAMIN = smallest possible value of ABS(SLAM)

    slamin = 0;

    for (i = 1; i <= fNpar; ++i) {

	if (step[i-1] != 0) {

           ratio = TMath_Abs(start[i-1] / step[i-1]);

           if (slamin == 0)    slamin = ratio;

           if (ratio < slamin) slamin = ratio;

        }

	fX[i-1] = start[i-1] + step[i-1];

    }

    if (slamin == 0) slamin = fEpsmac;

    slamin *= fEpsma2;

    nparx = fNpar;



    mninex(fX);

    (*fFCN)(nparx, fGin, f1, fU, 4);    ++fNfcn;

    ++nxypt;

    chpq[nxypt-1] = charal[nxypt-1];

    xpq[nxypt-1] = 1;

    ypq[nxypt-1] = f1;

    if (f1 < fstart) {

	fvmin = f1;

	xvmin = 1;

    }

//*-*-                        . quadr interp using slope GDEL and two points

    slam   = 1;

    toler8 = toler;

    slamax = 5;

    flast  = f1;

//*-*-                        can iterate on two-points (cut) if no imprvmnt



    do {

       denom = (flast - fstart - slope*slam)*2 / (slam*slam);

       slam  = 1;

       if (denom != 0)    slam = -slope / denom;

       if (slam < 0)      slam = slamax;

       if (slam > slamax) slam = slamax;

       if (slam < toler8) slam = toler8;

       if (slam < slamin) {

          l80 = 1;

          break;

       }

       if (TMath_Abs(slam - 1) < toler8 && f1 < fstart) {

          l70 = 1;

          break;

       }

       if (TMath_Abs(slam - 1) < toler8) slam = toler8 + 1;

       if (nxypt >= 12) {

         l65 = 1;

         break;

       }

       for (i = 1; i <= fNpar; ++i) { fX[i-1] = start[i-1] + slam*step[i-1]; }

       mninex(fX);

       (*fFCN)(fNpar, fGin, f2, fU, 4);    ++fNfcn;

       ++nxypt;

       chpq[nxypt-1] = charal[nxypt-1];

       xpq[nxypt-1]  = slam;

       ypq[nxypt-1]  = f2;

       if (f2 < fvmin) {

          fvmin = f2;

          xvmin = slam;

       }

       if (fstart == fvmin) {

          flast  = f2;

          toler8 = toler*slam;

          overal = slam - toler8;

          slamax = overal;

       }

    } while (fstart == fvmin);



    if (!l65 && !l70 && !l80) {

//*-*-                                       . quadr interp using 3 points

       xvals[0] = xpq[0];

       fvals[0] = ypq[0];

       xvals[1] = xpq[nxypt-2];

       fvals[1] = ypq[nxypt-2];

       xvals[2] = xpq[nxypt-1];

       fvals[2] = ypq[nxypt-1];

//*-*-                            begin iteration, calculate desired step

       do {

          slamax = TMath_Max(slamax,TMath_Abs(xvmin)*2);

          mnpfit(xvals, fvals, 3, coeff, sdev);

          if (coeff[2] <= 0) {

             slopem = coeff[2]*2*xvmin + coeff[1];

             if (slopem <= 0)  slam = xvmin + slamax;

             else              slam = xvmin - slamax;

          } else {

             slam = -coeff[1] / (coeff[2]*2);

             if (slam > xvmin + slamax) slam = xvmin + slamax;

             if (slam < xvmin - slamax) slam = xvmin - slamax;

          }

          if (slam > 0) if (slam > overal) slam = overal;

          else          if (slam < undral) slam = undral;



//*-*-              come here if step was cut below

          do {

             toler9 = TMath_Max(toler8,TMath_Abs(toler8*slam));

             for (ipt = 1; ipt <= 3; ++ipt) {

                if (TMath_Abs(slam - xvals[ipt-1]) < toler9) {

                   l70 = 1;

                   break;

                }

             }

             if (l70) break;

//*-*-               take the step

             if (nxypt >= 12) {

                l65 = 1;

                break;

             }

             for (i = 1; i <= fNpar; ++i) { fX[i-1] = start[i-1] + slam*step[i-1]; }

             mninex(fX);

             (*fFCN)(nparx, fGin, f3, fU, 4);    ++fNfcn;

             ++nxypt;

             chpq[nxypt-1] = charal[nxypt-1];

             xpq[nxypt-1]  = slam;

             ypq[nxypt-1]  = f3;

//*-*-            find worst previous point out of three

             fvmax = fvals[0];

             nvmax = 1;

             if (fvals[1] > fvmax) {

                fvmax = fvals[1];

                nvmax = 2;

             }

             if (fvals[2] > fvmax) {

                fvmax = fvals[2];

                nvmax = 3;

             }

//*-*-             if latest point worse than all three previous, cut step

             if (f3 >= fvmax) {

                if (nxypt >= 12) {

                   l65 = 1;

                   break;

                }

                if (slam > xvmin) overal = TMath_Min(overal,slam - toler8);

                if (slam < xvmin) undral = TMath_Max(undral,slam + toler8);

                slam = (slam + xvmin)*.5;

             }

          } while (f3 >= fvmax);



//*-*-             prepare another iteration, replace worst previous point

          if (l65 || l70) break;



          xvals[nvmax-1] = slam;

          fvals[nvmax-1] = f3;

          if (f3 < fvmin) {

             fvmin = f3;

             xvmin = slam;

          } else {

             if (slam > xvmin) overal = TMath_Min(overal,slam - toler8);

             if (slam < xvmin) undral = TMath_Max(undral,slam + toler8);

          }

       } while (nxypt < 12);

    }



//*-*-                                           . . end of iteration . . .

//*-*-           stop because too many iterations

    if (!l70 && !l80) {

       cmess = " LINE SEARCH HAS EXHAUSTED THE LIMIT OF FUNCTION CALLS ";

       if (ldebug) {

          fPrintf(" MNLINE DEBUG: steps=");

          for (kk = 1; kk <= fNpar; ++kk) {

             fPrintf("  %12.4g",step[kk-1]);

          }

       }

    }

//*-*-           stop because within tolerance

    if (l70) cmess = " LINE SEARCH HAS ATTAINED TOLERANCE ";

    if (l80) cmess = " STEP SIZE AT ARITHMETICALLY ALLOWED MINIMUM";



    fAmin = fvmin;

    for (i = 1; i <= fNpar; ++i) {

	fDirin[i-1] = step[i-1]*xvmin;

	fX[i-1]     = start[i-1] + fDirin[i-1];

    }

    mninex(fX);

    if (xvmin < 0) {

	mnwarn("D", "MNLINE", " LINE MINIMUM IN BACKWARDS DIRECTION");

    }

    if (fvmin == fstart) {

	mnwarn("D", "MNLINE", " LINE SEARCH FINDS NO IMPROVEMENT ");

    }

    if (ldebug) {

        fPrintf(" AFTER%3d POINTS,%s",nxypt,(const char*)cmess);

	mnplot(xpq, ypq, chpq, nxypt, fNpagwd, fNpagln);

    }

} /* mnline_ */



//______________________________________________________________________________

void Midnight::mnmatu(MInt kode)

{

//*-*-*-*-*-*-*-*Prints the covariance matrix v when KODE=1*-*-*-*-*-*-*-*-*

//*-*            ==========================================

//*-*        always prints the global correlations, and

//*-*        calculates and prints the individual correlation coefficients

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    /* Local variables */

    static MDouble vline[kMAXDIM];

    static MInt ndex, i, j, m, n, ncoef, nparm, id, it, ix;

    static MInt nsofar, ndi, ndj, iso, isw2, isw5;

    static MString ctemp;



    isw2 = fISW[1];

    if (isw2 < 1) {

	fPrintf("%s",(const char*)fCovmes[isw2]);

	return;

    }

    if (fNpar == 0) {

	fPrintf(" MNMATU: NPAR=0");

	return;

    }

//*-*-                                      . . . . .external error matrix

    if (kode == 1) {

	isw5    = fISW[4];

	fISW[4] = 2;

	mnemat(fP, fMaxint);

	if (isw2 < 3) {

	    fPrintf("%s",(const char*)fCovmes[isw2]);

	}

	fISW[4] = isw5;

    }

//*-*-                                      . . . . . correlation coeffs. .

    if (fNpar <= 1) return;

    mnwerr();

//*-*-    NCOEF is number of coeff. that fit on one line, not to exceed 20

    ncoef = (fNpagwd - 19) / 6;

    ncoef = TMath_Min(ncoef,20);

    nparm = TMath_Min(fNpar,ncoef);

    fPrintf(" PARAMETER  CORRELATION COEFFICIENTS  ");

    ctemp = "       NO.  GLOBAL";

    for (id = 1; id <= nparm; ++id) {

        ctemp += Form(" %6d",fNexofi[id-1]);

    }

    fPrintf("%s",(const char*)ctemp);

    for (i = 1; i <= fNpar; ++i) {

	ix  = fNexofi[i-1];

	ndi = i*(i + 1) / 2;

	for (j = 1; j <= fNpar; ++j) {

	    m    = TMath_Max(i,j);

	    n    = TMath_Min(i,j);

	    ndex = m*(m-1) / 2 + n;

	    ndj  = j*(j + 1) / 2;

	    vline[j-1] = fVhmat[ndex-1] / TMath_Sqrt(TMath_Abs(fVhmat[ndi-1]*fVhmat[ndj-1]));

	}

	nparm = TMath_Min(fNpar,ncoef);

        ctemp = Form("      %3d  %7.5f ",ix);

	for (it = 1; it <= nparm; ++it) {

            ctemp += Form(" %6.3f",vline[it-1]);

	}

        fPrintf("%s",(const char*)ctemp);

	if (i <= nparm) continue;

        ctemp = "                   ";

	for (iso = 1; iso <= 10; ++iso) {

	    nsofar = nparm;

	    nparm  = TMath_Min(fNpar,nsofar + ncoef);

	    for (it = nsofar + 1; it <= nparm; ++it) {

                ctemp += Form(" %6.3f",vline[it-1]);

	    }

            fPrintf("%s",(const char*)ctemp);

	    if (i <= nparm) break;

	}

    }

    if (isw2 < 3) {

        fPrintf(" %s",(const char*)fCovmes[isw2]);

    }

} /* mnmatu_ */



//______________________________________________________________________________

void Midnight::mnmigr()

{

//*-*-*-*-*-*-*-*-*Performs a local function minimization*-*-*-*-*-*-*-*-*-*

//*-*              ======================================

//*-*        Performs a local function minimization using basically the

//*-*        method of Davidon-Fletcher-Powell as modified by Fletcher

//*-*        ref. -- Fletcher, Comp.J. 13,317 (1970)   "switching method"

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    /* Local variables */

    static MDouble gdel, gami, flnu[kMAXDIM], vlen, step[kMAXDIM], dsum, gssq, vsum, d;

    static MDouble fzero, gs[kMAXDIM], fs, ri, vg[kMAXDIM], delgam, rhotol;

    static MDouble gdgssq, gvg, vgi, xxs[kMAXDIM];

    static MInt npfn, ndex, iext, i, j, m, n, npsdf, nparx;

    static MInt iswtr, lined2, kk, nfcnmg, nrstrt,iter;

    static MBool ldebug;

    static MDouble toler = 0.05;



    if (fNpar <= 0) return;

    if (fAmin == fUndefi) mnamin();

    ldebug  = kFALSE; if ( fIdbg[4] >= 1) ldebug = kTRUE;

    fCfrom  = "MIGRAD  ";

    fNfcnfr = fNfcn;

    nfcnmg  = fNfcn;

    fCstatu = "INITIATE  ";

    iswtr   = fISW[4] - 2*fItaur;

    npfn    = fNfcn;

    nparx   = fNpar;

    vlen    = (MDouble) (fNpar*(fNpar + 1) / 2);

    nrstrt  = 0;

    npsdf   = 0;

    lined2  = 0;

    fISW[3] = -1;

    rhotol  = fApsi*.001;

    if (iswtr >= 1) {

	fPrintf(" START MIGRAD MINIMIZATION.  STRATEGY%2d.  CONVERGENCE WHEN EDM .LT.%9.2e",fIstrat,rhotol);

    }

//*-*-                                          initialization strategy

    if (fIstrat < 2 || fISW[1] >= 3) goto L2;

//*-*-                               come (back) here to restart completely

L1:

    if (nrstrt > fIstrat) {

	fCstatu = "FAILED    ";

	fISW[3] = -1;

	goto L230;

    }

//*-*-                                     . get full covariance and gradient

    mnhess();

    mnwerr();

    npsdf = 0;

    if (fISW[1] >= 1) goto L10;

//*-*-                                       . get gradient at start point

L2:

    mninex(fX);

    if (fISW[2] == 1) {

	(*fFCN)(nparx, fGin, fzero, fU, 2);	++fNfcn;

    }

    mnderi();

    if (fISW[1] >= 1) goto L10;

//*-*-                                  sometimes start with diagonal matrix

    for (i = 1; i <= fNpar; ++i) {

	xxs[i-1]  = fX[i-1];

	step[i-1] = 0;

    }

//*-*-                          do line search if second derivative negative

    ++lined2;

    if (lined2 < (fIstrat + 1)*fNpar) {

	for (i = 1; i <= fNpar; ++i) {

	    if (fG2[i-1] > 0) continue;

	    if (fGrd[i-1] > 0) step[i-1] = -TMath_Abs(fGstep[i-1]);

	    else               step[i-1] =  TMath_Abs(fGstep[i-1]);

	    gdel = step[i-1]*fGrd[i-1];

	    fs   = fAmin;

	    mnline(xxs, fs, step, gdel, toler);

	    mnwarn("D", "MNMIGR", "Negative G2 line search");

	    iext = fNexofi[i-1];

	    if (ldebug) {

		fPrintf(" Negative G2 line search, param %3d %13.3g%13.3g",iext,fs,fAmin);

	    }

	    goto L2;

	}

    }

//*-*-                          make diagonal error matrix

    for (i = 1; i <= fNpar; ++i) {

	ndex = i*(i-1) / 2;

	for (j = 1; j <= i-1; ++j) {

	    ++ndex;

	    fVhmat[ndex-1] = 0;

	}

	++ndex;

	if (fG2[i-1] <= 0) fG2[i-1] = 1;

	fVhmat[ndex-1] = 2 / fG2[i-1];

    }

    fDcovar = 1;

    if (ldebug) {

	fPrintf(" DEBUG MNMIGR, STARTING MATRIX DIAGONAL,  VHMAT=");

	for (kk = 1; kk <= MInt(vlen); ++kk) {

	    fPrintf(" %10.2g",fVhmat[kk-1]);

	}

    }

//*-*-                                        ready to start first iteration

L10:

    ++nrstrt;

    if (nrstrt > fIstrat + 1) {

	fCstatu = "FAILED    ";

	goto L230;

    }

    fs = fAmin;

//*-*-                                       . . . get EDM and set up loop

    fEDM = 0;

    for (i = 1; i <= fNpar; ++i) {

	gs[i-1]  = fGrd[i-1];

	xxs[i-1] = fX[i-1];

	ndex     = i*(i-1) / 2;

	for (j = 1; j <= i-1; ++j) {

	    ++ndex;

	    fEDM += gs[i-1]*fVhmat[ndex-1]*gs[j-1];

	}

	++ndex;

	fEDM += gs[i-1]*gs[i-1]*.5*fVhmat[ndex-1];

    }

    fEDM = fEDM*.5*(fDcovar*3 + 1);

    if (fEDM < 0) {

	mnwarn("W", "MIGRAD", "STARTING MATRIX NOT POS-DEFINITE.");

	fISW[1] = 0;

	fDcovar = 1;

	goto L2;

    }

    if (fISW[1] == 0) fEDM = fBigedm;

    iter = 0;

    mninex(fX);

    mnwerr();

    if (iswtr >= 1) mnprin(3, fAmin);

    if (iswtr >= 2) mnmatu(0);

//*-*-                                       . . . . .  start main loop

L24:

    if (fNfcn - npfn >= fNfcnmx) goto L190;

    gdel = 0;

    gssq = 0;

    for (i = 1; i <= fNpar; ++i) {

	ri = 0;

	gssq += gs[i-1]*gs[i-1];

	for (j = 1; j <= fNpar; ++j) {

	    m    = TMath_Max(i,j);

	    n    = TMath_Min(i,j);

	    ndex = m*(m-1) / 2 + n;

	    ri  += fVhmat[ndex-1]*gs[j-1];

	}

	step[i-1] = ri*-.5;

	gdel += step[i-1]*gs[i-1];

    }

    if (gssq == 0) {

	mnwarn("D", "MIGRAD", " FIRST DERIVATIVES OF FCN ARE ALL ZERO");

	goto L300;

    }

//*-*-                if gdel positive, V not posdef

    if (gdel >= 0) {

	mnwarn("D", "MIGRAD", " NEWTON STEP NOT DESCENT.");

	if (npsdf == 1) goto L1;

	mnpsdf();

	npsdf = 1;

	goto L24;

    }

//*-*-                                       . . . . do line search

    mnline(xxs, fs, step, gdel, toler);

    if (fAmin == fs) goto L200;

    fCfrom  = "MIGRAD  ";

    fNfcnfr = nfcnmg;

    fCstatu = "PROGRESS  ";

//*-*-                                       . get gradient at new point

    mninex(fX);

    if (fISW[2] == 1) {

	(*fFCN)(nparx, fGin, fzero, fU, 2);	++fNfcn;

    }

    mnderi();

//*-*-                                        . calculate new EDM

    npsdf = 0;

L81:

    fEDM = 0;

    gvg = 0;

    delgam = 0;

    gdgssq = 0;

    for (i = 1; i <= fNpar; ++i) {

	ri  = 0;

	vgi = 0;

	for (j = 1; j <= fNpar; ++j) {

	    m    = TMath_Max(i,j);

	    n    = TMath_Min(i,j);

	    ndex = m*(m-1) / 2 + n;

	    vgi += fVhmat[ndex-1]*(fGrd[j-1] - gs[j-1]);

	    ri  += fVhmat[ndex-1]*fGrd[j-1];

	}

	vg[i-1] = vgi*.5;

	gami    = fGrd[i-1] - gs[i-1];

	gdgssq += gami*gami;

	gvg    += gami*vg[i-1];

	delgam += fDirin[i-1]*gami;

	fEDM   += fGrd[i-1]*ri*.5;

    }

    fEDM = fEDM*.5*(fDcovar*3 + 1);

//*-*-                         . if EDM negative,  not positive-definite

    if (fEDM < 0 || gvg <= 0) {

	mnwarn("D", "MIGRAD", "NOT POS-DEF. EDM OR GVG NEGATIVE.");

	fCstatu = "NOT POSDEF";

	if (npsdf == 1) goto L230;

	mnpsdf();

	npsdf = 1;

	goto L81;

    }

//*-*-                           print information about this iteration

    ++iter;

    if (iswtr >= 3 || iswtr == 2 && iter % 10 == 1) {

	mnwerr();

	mnprin(3, fAmin);

    }

    if (gdgssq == 0) {

	mnwarn("D", "MIGRAD", "NO CHANGE IN FIRST DERIVATIVES OVER LAST STEP");

    }

    if (delgam < 0) {

	mnwarn("D", "MIGRAD", "FIRST DERIVATIVES INCREASING ALONG SEARCH LINE");

    }

//*-*-                                       .  update covariance matrix

    fCstatu = "IMPROVEMNT";

    if (ldebug) {

	fPrintf(" VHMAT 1 =");

	for (kk = 1; kk <= 10; ++kk) {

            fPrintf(" %10.2g",fVhmat[kk-1]);

	}

    }

    dsum = 0;

    vsum = 0;

    for (i = 1; i <= fNpar; ++i) {

	for (j = 1; j <= i; ++j) {

	    d = fDirin[i-1]*fDirin[j-1] / delgam - vg[i-1]*vg[j-1] / gvg;

	    dsum += TMath_Abs(d);

	    ndex  = i*(i-1) / 2 + j;

	    fVhmat[ndex-1] += d*2;

	    vsum += TMath_Abs(fVhmat[ndex-1]);

	}

    }

//*-*-               smooth local fluctuations by averaging DCOVAR

    fDcovar = (fDcovar + dsum / vsum)*.5;

    if (iswtr >= 3 || ldebug) {

	fPrintf(" RELATIVE CHANGE IN COV. MATRIX=%5.1f per cent",fDcovar*100);

    }

    if (ldebug) {

	fPrintf(" VHMAT 2 =");

	for (kk = 1; kk <= 10; ++kk) {

            fPrintf(" %10.3g",fVhmat[kk-1]);

	}

    }

    if (delgam <= gvg) goto L135;

    for (i = 1; i <= fNpar; ++i) {

	flnu[i-1] = fDirin[i-1] / delgam - vg[i-1] / gvg;

    }

    for (i = 1; i <= fNpar; ++i) {

	for (j = 1; j <= i; ++j) {

	    ndex = i*(i-1) / 2 + j;

	    fVhmat[ndex-1] += gvg*2*flnu[i-1]*flnu[j-1];

	}

    }

L135:

//*-*-                                             and see if converged

    if (fEDM < rhotol*.1) goto L300;

//*-*-                                   if not, prepare next iteration

    for (i = 1; i <= fNpar; ++i) {

	xxs[i-1] = fX[i-1];

	gs[i-1]  = fGrd[i-1];

    }

    fs = fAmin;

    if (fISW[1] == 0 && fDcovar < .5)  fISW[1] = 1;

    if (fISW[1] == 3 && fDcovar > .1)  fISW[1] = 1;

    if (fISW[1] == 1 && fDcovar < .05) fISW[1] = 3;

    goto L24;

//*-*-                                       . . . . .  end main loop

//*-*-                                        . . call limit in MNMIGR

L190:

    fISW[0] = 1;

    if (fISW[4] >= 0) {

	fPrintf(" CALL LIMIT EXCEEDED IN MIGRAD.");

    }

    fCstatu = "CALL LIMIT";

    goto L230;

//*-*-                                        . . fails to improve . .

L200:

    if (iswtr >= 1) {

	fPrintf(" MIGRAD FAILS TO FIND IMPROVEMENT");

    }

    for (i = 1; i <= fNpar; ++i) { fX[i-1] = xxs[i-1]; }

    if (fEDM < rhotol) goto L300;

    if (fEDM < TMath_Abs(fEpsma2*fAmin)) {

	if (iswtr >= 0) {

	    fPrintf(" MACHINE ACCURACY LIMITS FURTHER IMPROVEMENT.");

	}

	goto L300;

    }

    if (fIstrat < 1) {

	if (fISW[4] >= 0) {

	    fPrintf(" MIGRAD FAILS WITH STRATEGY=0.   WILL TRY WITH STRATEGY=1.");

	}

	fIstrat = 1;

    }

    goto L1;

//*-*-                                        . . fails to converge

L230:

    if (iswtr >= 0) {

	fPrintf(" MIGRAD TERMINATED WITHOUT CONVERGENCE.");

    }

    if (fISW[1] == 3) fISW[1] = 1;

    fISW[3] = -1;

    goto L400;

//*-*-                                        . . apparent convergence

L300:

    if (iswtr >= 0) {

	fPrintf(" MIGRAD MINIMIZATION HAS CONVERGED.");

    }

    if (fItaur == 0) {

	if (fIstrat >= 2 || (fIstrat == 1 && fISW[1] < 3)) {

	    if (fISW[4] >= 0) {

		fPrintf(" MIGRAD WILL VERIFY CONVERGENCE AND ERROR MATRIX.");

	    }

	    mnhess();

	    mnwerr();

	    npsdf = 0;

	    if (fEDM > rhotol) goto L10;

	}

    }

    fCstatu = "CONVERGED ";

    fISW[3] = 1;

//*-*-                                          come here in any case

L400:

    fCfrom  = "MIGRAD  ";

    fNfcnfr = nfcnmg;

    mninex(fX);

    mnwerr();

    if (iswtr >= 0) mnprin(3, fAmin);

    if (iswtr >= 1) mnmatu(1);

} /* mnmigr_ */



//______________________________________________________________________________

void Midnight::mnmnos()

{

//*-*-*-*-*-*-*-*-*-*-*Performs a MINOS error analysis*-*-*-*-*-*-*-*-*-*-*-*

//*-*                  ===============================

//*-*        Performs a MINOS error analysis on those parameters for

//*-*        which it is requested on the MINOS command by calling

//*-*        MNMNOT for each parameter requested.

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    /* Local variables */

    static MDouble val2mi, val2pl;

    static MInt nbad, ilax, ilax2, ngood, nfcnmi, iin, knt;



    if (fNpar <= 0) goto L700;

    ngood = 0;

    nbad = 0;

    nfcnmi = fNfcn;

//*-*-                                     . loop over parameters requested

    for (knt = 1; knt <= fNpar; ++knt) {

	if (MInt(fWord7[1]) == 0) {

	    ilax = fNexofi[knt-1];

	} else {

	    if (knt >= 7) break;

	    ilax = MInt(fWord7[knt]);

	    if (ilax == 0) break;

	    if (ilax > 0 && ilax <= fNu) {

		if (fNiofex[ilax-1] > 0) goto L565;

	    }

            fPrintf(" PARAMETER NUMBER %3d NOT VARIABLE. IGNORED.",ilax);

	    continue;

	}

L565:

//*-*-                                        calculate one pair of M E s

	ilax2 = 0;

	mnmnot(ilax, ilax2, val2pl, val2mi);

	if (fLnewmn) goto L650;

//*-*-                                         update NGOOD and NBAD

	iin = fNiofex[ilax-1];

	if (fErp[iin-1] > 0) ++ngood;

	else                   ++nbad;

	if (fErn[iin-1] < 0) ++ngood;

	else                   ++nbad;

    }

//*-*-                                          end of loop . . . . . . .

//*-*-                                       . . . . printout final values .

    fCfrom  = "MINOS   ";

    fNfcnfr = nfcnmi;

    fCstatu = "UNCHANGED ";

    if (ngood == 0 && nbad == 0) goto L700;

    if (ngood > 0 && nbad == 0)  fCstatu = "SUCCESSFUL";

    if (ngood == 0 && nbad > 0)  fCstatu = "FAILURE   ";

    if (ngood > 0 && nbad > 0)   fCstatu = "PROBLEMS  ";

    if (fISW[4] >= 0)    mnprin(4, fAmin);

    if (fISW[4] >= 2)    mnmatu(0);

    return;

//*-*-                                       . . . new minimum found . . . .

L650:

    fCfrom  = "MINOS   ";

    fNfcnfr = nfcnmi;

    fCstatu = "NEW MINIMU";

    if (fISW[4] >= 0) mnprin(4, fAmin);

    fPrintf(" NEW MINIMUM FOUND.  GO BACK TO MINIMIZATION STEP.");

    fPrintf(" =================================================");

    fPrintf("                                                  V");

    fPrintf("                                                  V");

    fPrintf("                                                  V");

    fPrintf("                                               VVVVVVV");

    fPrintf("                                                VVVVV");

    fPrintf("                                                 VVV");

    fPrintf("                                                  V");

    fPrintf("");

    return;

L700:

    fPrintf(" THERE ARE NO MINOS ERRORS TO CALCULATE.");

} /* mnmnos_ */



//______________________________________________________________________________

void Midnight::mnmnot(MInt ilax, MInt ilax2, MDouble &val2pl, MDouble &val2mi)

{

//*-*-*-*-*-*Performs a MINOS error analysis on one parameter*-*-*-*-*-*-*-*-*

//*-*        ================================================

//*-*        The parameter ILAX is varied, and the minimum of the

//*-*        function with respect to the other parameters is followed

//*-*        until it crosses the value FMIN+UP.

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    /* System generated locals */

    MInt i__1;



    /* Local variables */

    static MDouble xdev[kMAXDIM], delu, aopt, eros;

    static MDouble w[kMAXDIM], abest, xunit, dc, ut, sigsav, du1;

    static MDouble fac, gcc[kMAXDIM], sig, sav;

    static MInt marc, isig, mpar, ndex, imax, indx, ierr, i, j;

    static MInt iercr, it, istrav, nfmxin, nlimit, isw2, isw4;

    static MString csig;



//*-*-                                       . . save and prepare start vals

    isw2    = fISW[1];

    isw4    = fISW[3];

    sigsav  = fEDM;

    istrav  = fIstrat;

    dc      = fDcovar;

    fLnewmn = kFALSE;

    fApsi   = fEpsi*.5;

    abest   = fAmin;

    mpar    = fNpar;

    nfmxin  = fNfcnmx;

    for (i = 1; i <= mpar; ++i) { fXt[i-1] = fX[i-1]; }

    i__1 = mpar*(mpar + 1) / 2;

    for (j = 1; j <= i__1; ++j) { fVthmat[j-1] = fVhmat[j-1]; }

    for (i = 1; i <= mpar; ++i) {

	gcc[i-1] = fGlobcc[i-1];

	w[i-1]   = fWerr[i-1];

    }

    it = fNiofex[ilax-1];

    fErp[it-1] = 0;

    fErn[it-1] = 0;

    mninex(fXt);

    ut = fU[ilax-1];

    if (fNvarl[ilax-1] == 1) {

	fAlim[ilax-1] = ut - w[it-1]*100;

	fBlim[ilax-1] = ut + w[it-1]*100;

    }

    ndex  = it*(it + 1) / 2;

    xunit = TMath_Sqrt(fUp / fVthmat[ndex-1]);

    marc  = 0;

    for (i = 1; i <= mpar; ++i) {

	if (i == it) continue;

	++marc;

	imax = TMath_Max(it,i);

	indx = imax*(imax-1) / 2 + TMath_Min(it,i);

	xdev[marc-1] = xunit*fVthmat[indx-1];

    }

//*-*-                          fix the parameter in question

    mnfixp(it-1, ierr);

    if (ierr > 0) {

	fPrintf(" MINUIT ERROR. CANNOT FIX PARAMETER%4d    INTERNAL%3d",ilax,it);

	goto L700;

    }

//*-*-                      . . . . . Nota Bene: from here on, NPAR=MPAR-1

//*-*-     Remember: MNFIXP squeezes IT out of X, XT, WERR, and VHMAT,

//*-*-                                                   not W, VTHMAT

    for (isig = 1; isig <= 2; ++isig) {

	if (isig == 1) {

	    sig  = 1;

	    csig = "POSI";

	} else {

	    sig  = -1;

	    csig = "NEGA";

	}

//*-*-                                       . sig=sign of error being calcd

	if (fISW[4] > 1) {

	    fPrintf(" DETERMINATION OF %sTIVE MINOS ERROR FOR PARAMETER%d"

                            ,(const char*)csig,ilax

                            ,(const char*)fCpnam[ilax-1]);

	}

	if (fISW[1] <= 0) {

	    mnwarn("D", "MINOS", "NO COVARIANCE MATRIX.");

	}

	nlimit     = fNfcn + nfmxin;

	fIstrat    = TMath_Max(istrav-1,0);

	du1        = w[it-1];

	fU[ilax-1] = ut + sig*du1;

	fU[ilax-1] = TMath_Min(fU[ilax-1],fBlim[ilax-1]);

	fU[ilax-1] = TMath_Max(fU[ilax-1],fAlim[ilax-1]);

	delu = fU[ilax-1] - ut;

//*-*-        stop if already at limit with negligible step size

	if (TMath_Abs(delu) / (TMath_Abs(ut) + TMath_Abs(fU[ilax-1])) < fEpsmac) goto L440;

	fac = delu / w[it-1];

	for (i = 1; i <= fNpar; ++i) {

	    fX[i-1] = fXt[i-1] + fac*xdev[i-1];

	}

	if (fISW[4] > 1) {

	    fPrintf(" PARAMETER%4d SET TO%11.3e + %10.3e = %12.3e",ilax,ut,delu,fU[ilax-1]);

	}

//*-*-                                       loop to hit AMIN+UP

	fKe1cr  = ilax;

	fKe2cr  = 0;

	fXmidcr = fU[ilax-1];

	fXdircr = delu;



	fAmin = abest;

	fNfcnmx = nlimit - fNfcn;

	mncros(aopt, iercr);

	if (abest - fAmin > fUp*.01) goto L650;

	if (iercr == 1) goto L440;

	if (iercr == 2) goto L450;

	if (iercr == 3) goto L460;

//*-*-                                       . error successfully calculated

	eros = fXmidcr - ut + aopt*fXdircr;

	if (fISW[4] > 1) {

	    fPrintf("         THE %4sTIVE MINOS ERROR OF PARAMETER%3d  %10s, IS %12.4e"

                           ,(const char*)csig,ilax

                           ,(const char*)fCpnam[ilax-1],eros);

	}

	goto L480;

//*-*-                                       . . . . . . . . failure returns

L440:

	if (fISW[4] >= 1) {

	    fPrintf("     THE %4sTIVE MINOS ERROR OF PARAMETER%3d, %s EXCEEDS ITS LIMIT."

                              ,(const char*)csig,ilax

                              ,(const char*)fCpnam[ilax-1]);

	}

	eros = fUndefi;

	goto L480;

L450:

	if (fISW[4] >= 1) {

	    fPrintf("         THE %4sTIVE MINOS ERROR%4d REQUIRES MORE THAN%5d FUNCTION CALLS."

                         ,(const char*)csig,ilax,nfmxin);

	}

	eros = 0;

	goto L480;

L460:

	if (fISW[4] >= 1) {

	    fPrintf("                         %4sTIVE MINOS ERROR NOT CALCULATED FOR PARAMETER%d"

                         ,(const char*)csig,ilax);

	}

	eros = 0;



L480:

	if (fISW[4] > 1) {

	    fPrintf("     **************************************************************************");

	}

	if (sig < 0) {

	    fErn[it-1] = eros;

	    if (ilax2 > 0 && ilax2 <= fNu) val2mi = fU[ilax2-1];

	} else {

	    fErp[it-1] = eros;

	    if (ilax2 > 0 && ilax2 <= fNu) val2pl = fU[ilax2-1];

	}

    }

//*-*-                                       . . parameter finished. reset v

//*-*-                      normal termination */

    fItaur = 1;

    mnfree(1);

    i__1 = mpar*(mpar + 1) / 2;

    for (j = 1; j <= i__1; ++j) { fVhmat[j-1] = fVthmat[j-1]; }

    for (i = 1; i <= mpar; ++i) {

	fWerr[i-1]   = w[i-1];

	fGlobcc[i-1] = gcc[i-1];

	fX[i-1]      = fXt[i-1];

    }

    mninex(fX);

    fEDM    = sigsav;

    fAmin   = abest;

    fISW[1] = isw2;

    fISW[3] = isw4;

    fDcovar = dc;

    goto L700;

//*-*-                      new minimum

L650:

    fLnewmn = kTRUE;

    fISW[1] = 0;

    fDcovar = 1;

    fISW[3] = 0;

    sav     = fU[ilax-1];

    fItaur  = 1;

    mnfree(1);

    fU[ilax-1] = sav;

    mnexin(fX);

    fEDM = fBigedm;

//*-*-                      in any case

L700:

    fItaur  = 0;

    fNfcnmx = nfmxin;

    fIstrat = istrav;

} /* mnmnot_ */



//______________________________________________________________________________

void Midnight::mnparm(MInt k1, MString cnamj, MDouble uk, MDouble wk, MDouble a, MDouble b, MInt &ierflg)

{

//*-*-*-*-*-*-*-*-*Implements one parameter definition*-*-*-*-*-*-*-*-*-*-*-*

//*-*              ===================================

//*-*        Called from MNPARS and user-callable

//*-*    Implements one parameter definition, that is:

//*-*          K     (external) parameter number

//*-*          CNAMK parameter name

//*-*          UK    starting value

//*-*          WK    starting step size or uncertainty

//*-*          A, B  lower and upper physical parameter limits

//*-*    and sets up (updates) the parameter lists.

//*-*    Output: IERFLG=0 if no problems

//*-*                  >0 if MNPARM unable to implement definition

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    /* Local variables */

    static MDouble vplu, a_small, gsmin, pinti, vminu, danger, sav, sav2;

    static MInt ierr, kint, in, ix, ktofix, lastin, kinfix, nvl;

    static MString cnamk, chbufi;



    MInt k = k1+1;

    cnamk   = cnamj;

    kint    = fNpar;

    if (k < 1 || k > fMaxext) {

//*-*-                    parameter number exceeds allowed maximum value

	fPrintf(" MINUIT USER ERROR.  PARAMETER NUMBER IS %3d  ALLOWED RANGE IS ONE TO %4d",k,fMaxext);

	goto L800;

    }

//*-*-                    normal parameter request

    ktofix = 0;

    if (fNvarl[k-1] < 0) goto L50;

//*-*-        previously defined parameter is being redefined

//*-*-                                    find if parameter was fixed

    for (ix = 1; ix <= fNpfix; ++ix) {

	if (fIpfix[ix-1] == k) ktofix = k;

    }

    if (ktofix > 0) {

	mnwarn("W", "PARAM DEF", "REDEFINING A FIXED PARAMETER.");

	if (kint >= fMaxint) {

	    fPrintf(" CANNOT RELEASE. MAX NPAR EXCEEDED.");

	    goto L800;

	}

	mnfree(1);

    }

//*-*-                      if redefining previously variable parameter

    if (fNiofex[k-1] > 0) kint = fNpar - 1;

L50:



//*-*-                                     . . .print heading

    if (fLphead && fISW[4] >= 0) {

	fPrintf(" PARAMETER DEFINITIONS:");

	fPrintf("    NO.   NAME         VALUE      STEP SIZE      LIMITS");

	fLphead = kFALSE;

    }

    if (wk > 0) goto L122;

//*-*-                                       . . .constant parameter . . . .

    if (fISW[4] >= 0) {

        fPrintf(" %5d %-10s %13.5e  constant",k,(const char*)cnamk,uk);

    }

    nvl = 0;

    goto L200;

L122:

    if (a == 0 && b == 0) {

//*-*-                                     variable parameter without limits

	nvl = 1;

	if (fISW[4] >= 0) {

            fPrintf(" %5d %-10s %13.5e%13.5e     no limits",k,(const char*)cnamk,uk,wk);

	}

    } else {

//*-*-                                        variable parameter with limits

	nvl = 4;

	fLnolim = kFALSE;

	if (fISW[4] >= 0) {

            fPrintf(" %5d '%-10s' %13.5e%13.5e  %13.5e%13.5e",k,(const char*)cnamk,uk,wk,a,b);

	}

    }

//*-*-                            . . request for another variable parameter

    ++kint;

    if (kint > fMaxint) {

	fPrintf(" MINUIT USER ERROR.   TOO MANY VARIABLE PARAMETERS.");

	goto L800;

    }

    if (nvl == 1) goto L200;

    if (a == b) {

	fPrintf(" USER ERROR IN MINUIT PARAMETER");

	fPrintf(" DEFINITION");

	fPrintf(" UPPER AND LOWER LIMITS EQUAL.");

	goto L800;

    }

    if (b < a) {

	sav = b;

	b = a;

	a = sav;

	mnwarn("W", "PARAM DEF", "PARAMETER LIMITS WERE REVERSED.");

	if (fLwarn) fLphead = kTRUE;

    }

    if (b - a > 1e7) {

	mnwarn("W", "PARAM DEF", Form("LIMITS ON PARAM%d TOO FAR APART.",k));

	if (fLwarn) fLphead = kTRUE;

    }

    danger = (b - uk)*(uk - a);

    if (danger < 0) {

	mnwarn("W", "PARAM DEF", "STARTING VALUE OUTSIDE LIMITS.");

    }

    if (danger == 0) {

	mnwarn("W", "PARAM DEF", "STARTING VALUE IS AT LIMIT.");

    }

L200:

//*-*-                          . . . input OK, set values, arrange lists,

//*-*-                                   calculate step sizes GSTEP, DIRIN

    fCfrom      = "PARAMETR";

    fNfcnfr     = fNfcn;

    fCstatu     = "NEW VALUES";

    fNu         = TMath_Max(fNu,k);

    fCpnam[k-1] = cnamk;

    fU[k-1]     = uk;

    fAlim[k-1]  = a;

    fBlim[k-1]  = b;

    fNvarl[k-1] = nvl;

    mnrset(1);

//*-*-                            K is external number of new parameter

//*-*-          LASTIN is the number of var. params with ext. param. no.< K

    lastin = 0;

    for (ix = 1; ix <= k-1; ++ix) { if (fNiofex[ix-1] > 0) ++lastin; }

//*-*-                KINT is new number of variable params, NPAR is old

    if (kint == fNpar) goto L280;

    if (kint > fNpar) {

//*-*-                         insert new variable parameter in list

	for (in = fNpar; in >= lastin + 1; --in) {

	    ix            = fNexofi[in-1];

	    fNiofex[ix-1] = in + 1;

	    fNexofi[in]   = ix;

	    fX[in]        = fX[in-1];

	    fXt[in]       = fXt[in-1];

	    fDirin[in]    = fDirin[in-1];

	    fG2[in]       = fG2[in-1];

	    fGstep[in]    = fGstep[in-1];

	}

    } else {

//*-*-                         remove variable parameter from list

	for (in = lastin + 1; in <= kint; ++in) {

	    ix            = fNexofi[in];

	    fNiofex[ix-1] = in;

	    fNexofi[in-1] = ix;

	    fX[in-1]      = fX[in];

	    fXt[in-1]     = fXt[in];

	    fDirin[in-1]  = fDirin[in];

	    fG2[in-1]     = fG2[in];

	    fGstep[in-1]  = fGstep[in];

	}

    }

L280:

    ix = k;

    fNiofex[ix-1] = 0;

    fNpar = kint;

//*-*-                                      lists are now arranged . . . .

    if (nvl > 0) {

	in            = lastin + 1;

	fNexofi[in-1] = ix;

	fNiofex[ix-1] = in;

	sav           = fU[ix-1];

	mnpint(sav, ix-1, pinti);

	fX[in-1]    = pinti;

	fXt[in-1]   = fX[in-1];

	fWerr[in-1] = wk;

	sav2        = sav + wk;

	mnpint(sav2, ix-1, pinti);

	vplu = pinti - fX[in-1];

	sav2 = sav - wk;

	mnpint(sav2, ix-1, pinti);

	vminu = pinti - fX[in-1];

	fDirin[in-1] = (TMath_Abs(vplu) + TMath_Abs(vminu))*.5;

	fG2[in-1] = fUp*2 / (fDirin[in-1]*fDirin[in-1]);

	gsmin = fEpsma2*8*TMath_Abs(fX[in-1]);

	fGstep[in-1] = TMath_Max(gsmin,fDirin[in-1]*.1);

	if (fAmin != fUndefi) {

	    a_small      = TMath_Sqrt(fEpsma2*(fAmin + fUp) / fUp);

	    fGstep[in-1] = TMath_Max(gsmin,a_small*fDirin[in-1]);

	}

	fGrd[in-1] = fG2[in-1]*fDirin[in-1];

//*-*-                  if parameter has limits

	if (fNvarl[k-1] > 1) {

	    if (fGstep[in-1] > .5) fGstep[in-1] = .5;

	    fGstep[in-1] = -fGstep[in-1];

	}

    }

    if (ktofix > 0) {

	kinfix = fNiofex[ktofix-1];

	if (kinfix > 0) mnfixp(kinfix-1, ierr);

	if (ierr > 0)   goto L800;

    }

    ierflg = 0;

    return;

//*-*-                  error on input, unable to implement request  . . . .

L800:

    ierflg = 1;

} /* mnparm_ */



//______________________________________________________________________________

void Midnight::mnpars(MString &crdbuf, MInt &icondn)

{

//*-*-*-*-*-*-*-*Implements one parameter definition*-*-*-*-*-*-*-*-*-*-*-*-*

//*-*            =========== =======================

//*-*        Called from MNREAD and user-callable

//*-*    Implements one parameter definition, that is:

//*-*       parses the string CRDBUF and calls MNPARM

//*-*

//*-* output conditions:

//*-*        ICONDN = 0    all OK

//*-*        ICONDN = 1    error, attempt to define parameter is ignored

//*-*        ICONDN = 2    end of parameter definitions

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    /* Local variables */

    static MDouble a, b, plist[30], fk, uk, wk, xk;

    static MInt ierr, kapo1, kapo2;

    static MInt k, llist, ibegin, lenbuf, istart, lnc, icy;

    static MString cnamk, comand, celmnt, ctemp;

    char stmp[128];



    lenbuf = strlen((const char*)crdbuf);

//*-*-                    find out whether fixed or free-field format

    kapo1 = strspn((const char*)crdbuf, "'");

    if (kapo1 == 0) goto L150;

    kapo2 = strspn((const char*)crdbuf + kapo1, "'");

    if (kapo2 == 0) goto L150;

//*-*-         new (free-field) format

    kapo2 += kapo1;

//*-*-                            skip leading blanks if any

    for (istart = 1; istart <= kapo1-1; ++istart) {

	if (crdbuf[istart-1] != ' ') goto L120;

    }

    goto L210;

L120:

//*-*-                              parameter number integer

    celmnt = crdbuf(istart-1, kapo1-istart);

    scanf((const char*)celmnt,fk);

    k = MInt(fk);

    if (k <= 0) goto L210;

    cnamk = "PARAM ";

    cnamk += celmnt;

    if (kapo2 - kapo1 > 1) {

	cnamk = crdbuf(kapo1, kapo2-1-kapo1);

    }

//*-*  special handling if comma or blanks and a comma follow 'name'

    for (icy = kapo2 + 1; icy <= lenbuf; ++icy) {

	if (crdbuf[icy-1] == ',') goto L139;

	if (crdbuf[icy-1] != ' ') goto L140;

    }

    uk = 0;

    wk = 0;

    a  = 0;

    b  = 0;

    goto L170;

L139:

    ++icy;

L140:

    ibegin = icy;

    ctemp = crdbuf(ibegin-1,lenbuf-ibegin);

    mncrck(ctemp, 20, comand, lnc, 30, plist, llist, ierr, fIsyswr);

    if (ierr > 0) goto L180;

    uk = plist[0];

    wk = 0;

    if (llist >= 2) wk = plist[1];

    a = 0;

    if (llist >= 3) a = plist[2];

    b = 0;

    if (llist >= 4) b = plist[3];

    goto L170;

//*-*-         old (fixed-field) format

L150:

    scanf((const char*)crdbuf,xk,stmp,uk,wk,a,b);

    cnamk = stmp;

    k = MInt(xk);

    if (k == 0)    goto L210;

//*-*-         parameter format cracked, implement parameter definition

L170:

    mnparm(k-1, cnamk, uk, wk, a, b, ierr);

    icondn = ierr;

    return;

//*-*-         format or other error

L180:

    icondn = 1;

    return;

//*-*-       end of data

L210:

    icondn = 2;

} /* mnpars_ */



//______________________________________________________________________________

void Midnight::mnpfit(MDouble *parx2p, MDouble *pary2p, MInt npar2p, MDouble *coef2p, MDouble &sdev2p)

{

//*-*-*-*-*-*-*-*-*-*To fit a parabola to npar2p points*-*-*-*-*-*-*-*-*-*-*

//*-*                ==================================

//*-*   npar2p   no. of points

//*-*   parx2p(i)   x value of point i

//*-*   pary2p(i)   y value of point i

//*-*

//*-*   coef2p(1...3)  coefficients of the fitted parabola

//*-*   y=coef2p(1) + coef2p(2)*x + coef2p(3)*x**2

//*-*   sdev2p= variance

//*-*   method : chi**2 = min equation solved explicitly

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    /* Local variables */

    static MDouble a, f, s, t, y, s2, x2, x3, x4, y2, cz[3], xm, xy, x2y;

    x2 = x3 = 0;

    MInt i;



    /* Parameter adjustments */

    --coef2p;

    --pary2p;

    --parx2p;



    /* Function Body */

    for (i = 1; i <= 3; ++i) { cz[i-1] = 0; }

    sdev2p = 0;

    if (npar2p < 3) goto L10;

    f = (MDouble) (npar2p);

//*-* --- center x values for reasons of machine precision

    xm  = 0;

    for (i = 1; i <= npar2p; ++i) { xm += parx2p[i]; }

    xm /= f;

    x2  = 0;

    x3  = 0;

    x4  = 0;

    y   = 0;

    y2  = 0;

    xy  = 0;

    x2y = 0;

    for (i = 1; i <= npar2p; ++i) {

	s    = parx2p[i] - xm;

	t    = pary2p[i];

	s2   = s*s;

	x2  += s2;

	x3  += s*s2;

	x4  += s2*s2;

	y   += t;

	y2  += t*t;

	xy  += s*t;

	x2y += s2*t;

    }

    a = (f*x4 - x2*x2)*x2 - f*(x3*x3);

    if (a == 0) goto L10;

    cz[2] = (x2*(f*x2y - x2*y) - f*x3*xy) / a;

    cz[1] = (xy - x3*cz[2]) / x2;

    cz[0] = (y - x2*cz[2]) / f;

    if (npar2p == 3) goto L6;

    sdev2p = y2 - (cz[0]*y + cz[1]*xy + cz[2]*x2y);

    if (sdev2p < 0) sdev2p = 0;

    sdev2p /= f - 3;

L6:

    cz[0] += xm*(xm*cz[2] - cz[1]);

    cz[1] -= xm*2*cz[2];

L10:

    for (i = 1; i <= 3; ++i) { coef2p[i] = cz[i-1]; }

} /* mnpfit_ */



//______________________________________________________________________________

void Midnight::mnpint(MDouble &pexti, MInt i1, MDouble &pinti)

{

//*-*-*-*-*-*-*Calculates the internal parameter value PINTI*-*-*-*-*-*-*-*

//*-*          =============================================

//*-*        corresponding  to the external value PEXTI for parameter I.

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    /* Local variables */

    static MDouble a, alimi, blimi, yy, yy2;

    static MInt igo;

    static MString chbuf2, chbufi;



    MInt i = i1+1;

    pinti   = pexti;

    igo     = fNvarl[i-1];

    if (igo == 4) {

//*-* --                          there are two limits

	alimi = fAlim[i-1];

	blimi = fBlim[i-1];

	yy = (pexti - alimi)*2 / (blimi - alimi) - 1;

	yy2 = yy*yy;

	if (yy2 >= 1 - fEpsma2) {

	    if (yy < 0) {

		a      = fVlimlo;

		chbuf2 = " IS AT ITS LOWER ALLOWED LIMIT.";

	    } else {

		a      = fVlimhi;

		chbuf2 = " IS AT ITS UPPER ALLOWED LIMIT.";

	    }

	    pinti   = a;

	    pexti   = alimi + (blimi - alimi)*.5*(TMath_Sin(a) + 1);

	    fLimset = kTRUE;

	    if (yy2 > 1) chbuf2 = " BROUGHT BACK INSIDE LIMITS.";

	    mnwarn("W", fCfrom, Form("VARIABLE%d%s",i,(const char*)chbuf2));

	} else {

	    pinti = TMath_ASin(yy);

	}

    }

} /* mnpint_ */



//______________________________________________________________________________

void Midnight::mnplot(MDouble *xpt, MDouble *ypt, MString *chpt, MInt nxypt, MInt npagwd, MInt npagln)

{

 //*-*-*-*Plots points in array xypt onto one page with labelled axes*-*-*-*-*

 //*-*    ===========================================================

 //*-*        NXYPT is the number of points to be plotted

 //*-*        XPT(I) = x-coord. of ith point

 //*-*        YPT(I) = y-coord. of ith point

 //*-*        CHPT(I) = character to be plotted at this position

 //*-*        the input point arrays XPT, YPT, CHPT are destroyed.

 //*-*

 //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



     /* Initialized data */



     static MString cdot   = ".";

     static MString cslash = "/";

     static MString cblank = " ";



     /* Local variables */

     static MDouble xmin, ymin, xmax, ymax, savx, savy, yprt;

     static MDouble bwidx, bwidy, xbest, ybest, ax, ay, bx, by;

     static MDouble xvalus[12], any, dxx, dyy;

     static MInt iten, i, j, k, maxnx, maxny, iquit, ni, linodd;

     static MInt nxbest, nybest, km1, ibk, isp1, nx, ny, ks, ix;

     static MString cline, chsav, chmess, chbest, ctemp;

     static MBool overpr;



     /* Function Body */

 //*-*  Computing MIN

     maxnx = TMath_Min(npagwd-20,100);

     if (maxnx < 10) maxnx = 10;

     maxny = npagln;

     if (maxny < 10) maxny = 10;

     if (nxypt <= 1) return;

     xbest  = xpt[0];

     ybest  = ypt[0];

     chbest = chpt[0];

 //*-*-        order the points by decreasing y

     km1 = nxypt - 1;

     for (i = 1; i <= km1; ++i) {

 	iquit = 0;

 	ni    = nxypt - i;

 	for (j = 1; j <= ni; ++j) {

 	    if (ypt[j-1] > ypt[j]) continue;

 	    savx     = xpt[j-1];

 	    xpt[j-1] = xpt[j];

 	    xpt[j]   = savx;

 	    savy     = ypt[j-1];

 	    ypt[j-1] = ypt[j];

 	    ypt[j]   = savy;

 	    chsav    = chpt[j-1];

 	    chpt[j-1]= chpt[j];

 	    chpt[j]  = chsav;

 	    iquit    = 1;

 	}

 	if (iquit == 0) break;

     }

 //*-*-        find extreme values

     xmax = xpt[0];

     xmin = xmax;

     for (i = 1; i <= nxypt; ++i) {

 	if (xpt[i-1] > xmax) xmax = xpt[i-1];

 	if (xpt[i-1] < xmin) xmin = xpt[i-1];

     }

     dxx   = (xmax - xmin)*.001;

     xmax += dxx;

     xmin -= dxx;

     mnbins(xmin, xmax, maxnx, xmin, xmax, nx, bwidx);

     ymax = ypt[0];

     ymin = ypt[nxypt-1];

     if (ymax == ymin) ymax = ymin + 1;

     dyy   = (ymax - ymin)*.001;

     ymax += dyy;

     ymin -= dyy;

     mnbins(ymin, ymax, maxny, ymin, ymax, ny, bwidy);

     any = (MDouble) ny;

 //*-*-        if first point is blank, it is an 'origin'

     if (chbest == cblank) goto L50;

     xbest = (xmax + xmin)*.5;

     ybest = (ymax + ymin)*.5;

 L50:

 //*-*-        find scale constants

     ax = 1 / bwidx;

     ay = 1 / bwidy;

     bx = -ax*xmin + 2;

     by = -ay*ymin - 2;

 //*-*-        convert points to grid positions

     for (i = 1; i <= nxypt; ++i) {

 	xpt[i-1] = ax*xpt[i-1] + bx;

 	ypt[i-1] = any - ay*ypt[i-1] - by;

     }

     nxbest = MInt((ax*xbest + bx));

     nybest = MInt((any - ay*ybest - by));

 //*-*-        print the points

     ny += 2;

     nx += 2;

     isp1 = 1;

     linodd = 1;

     overpr = kFALSE;



     for (i = 1; i <= ny; ++i) {

      cline.resize(nx+2);  // added Z.Sz.

 	for (ibk = 1; ibk <= nx; ++ibk) { cline[ibk-1] = ' '; }

// 	cline(nx) = '\0';  // removed Z.Sz.

 	cline(nx+1) = '\0';

 	cline(0)        = '.';

 	cline(nx-1)     = '.';

 	cline(nxbest-1) = '.';

 	if (i != 1 && i != nybest && i != ny) goto L320;

 	for (j = 1; j <= nx; ++j) { cline(j-1) = '.'; }

 L320:

 	yprt = ymax - MDouble(i-1)*bwidy;

 	if (isp1 > nxypt) goto L350;

 //*-*-        find the points to be plotted on this line

 	for (k = isp1; k <= nxypt; ++k) {

 	    ks = MInt(ypt[k-1]);

 	    if (ks > i) goto L345;

 	    ix = MInt(xpt[k-1]);

 	    if (cline(ix-1) == '.')   goto L340;

 	    if (cline(ix-1) == ' ') goto L340;

 	    if (cline(ix-1) == chpt[k-1](0))   continue;

 	    overpr = kTRUE;

 //*-*-        OVERPR is true if one or more positions contains more than

 //*-*-           one point

 	    cline(ix-1) = '&';

 	    continue;

 L340:

 	    cline(ix-1) = chpt[k-1](0);

 	}

 	isp1 = nxypt + 1;

 	goto L350;

 L345:

 	isp1 = k;

 L350:

 	if (linodd == 1 || i == ny) goto L380;

 	linodd = 1;

 	ctemp  = cline(0,nx);

 	fPrintf("                  %s",(const char*)ctemp);

 	goto L400;

 L380:

 	ctemp = cline(0,nx);

 	fPrintf(" %14.7g ..%s",yprt,(const char*)ctemp);

 	linodd = 0;

 L400:

 	;

     }

 //*-*-        print labels on x-axis every ten columns

     for (ibk = 1; ibk <= nx; ++ibk) {

 	cline[ibk-1] = ' ';

 	if (ibk % 10 == 1) cline[ibk-1] = '/';

     }

     fPrintf("                  %s",(const char*)cline);



     for (ibk = 1; ibk <= 12; ++ibk) {

 	xvalus[ibk-1] = xmin + MDouble(ibk-1)*10*bwidx;

     }

     fPrintf("           ");

     iten = (nx + 9) / 10;



// wrong: numbers should be printed in one line !!!!!!!!!!!!!!!!!!!!!!!!! Z.Sz.

     for (ibk = 1; ibk <= iten; ++ibk) {

         fPrintf(" %9.4g", xvalus[ibk-1]);

     }

// until here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



     chmess = " ";

     if (overpr) chmess = "   Overprint character is &";

     fPrintf("                         ONE COLUMN=%13.7g%s",bwidx,(const char*)chmess);

} /* mnplot_ */



//______________________________________________________________________________

void Midnight::mnpout(MInt iuext1, MString& chnam, MDouble &val, MDouble &err, MDouble &xlolim, MDouble &xuplim, MInt &iuint)

{

//*-*-*-*Provides the user with information concerning the current status*-*-*

//*-*    ================================================================

//*-*          of parameter number IUEXT. Namely, it returns:

//*-*        CHNAM: the name of the parameter

//*-*        VAL: the current (external) value of the parameter

//*-*        ERR: the current estimate of the parameter uncertainty

//*-*        XLOLIM: the lower bound (or zero if no limits)

//*-*        XUPLIM: the upper bound (or zero if no limits)

//*-*        IUINT: the internal parameter number (or zero if not variable,

//*-*           or negative if undefined).

//*-*  Note also:  If IUEXT is negative, then it is -internal parameter

//*-*           number, and IUINT is returned as the EXTERNAL number.

//*-*     Except for IUINT, this is exactly the inverse of MNPARM

//*-*     User-called

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    /* Local variables */

    static MInt iint, iext, nvl;



    MInt iuext = iuext1 + 1;

    xlolim = 0;

    xuplim = 0;

    err    = 0;

    if (iuext == 0) goto L100;

    if (iuext < 0) {

//*-*-                  internal parameter number specified

	iint  = -(iuext);

	if (iint > fNpar) goto L100;

	iext  = fNexofi[iint-1];

	iuint = iext;

    } else {

//*-*-                   external parameter number specified

	iext = iuext;

	if (iext == 0) goto L100;

	if (iext > fNu) goto L100;

	iint  = fNiofex[iext-1];

	iuint = iint;

    }

//*-*-                    in both cases

    nvl = fNvarl[iext-1];

    if (nvl < 0) goto L100;

    chnam = fCpnam[iext-1];

    val   = fU[iext-1];

    if (iint > 0) err = fWerr[iint-1];

    if (nvl == 4) {

	xlolim = fAlim[iext-1];

	xuplim = fBlim[iext-1];

    }

    return;

//*-*-               parameter is undefined

L100:

    iuint = -1;

    chnam = "undefined";

    val = 0;

} /* mnpout_ */



//______________________________________________________________________________

void Midnight::mnprin(MInt inkode, MDouble fval)

{

//*-*-*-*Prints the values of the parameters at the time of the call*-*-*-*-*

//*-*    ===========================================================

//*-*        also prints other relevant information such as function value,

//*-*        estimated distance to minimum, parameter errors, step sizes.

//*-*

//*-*         According to the value of IKODE, the printout is:/

//*-*    IKODE=INKODE= 0    only info about function value

//*-*                  1    parameter values, errors, limits

//*-*                  2    values, errors, step sizes, internal values

//*-*                  3    values, errors, step sizes, first derivs.

//*-*                  4    values, parabolic errors, MINOS errors

//*-*    when INKODE=5, MNPRIN chooses IKODE=1,2, or 3, according to ISW(2)

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    /* Initialized data */



    static MString cblank = "           ";

    static MString cnambf = "           ";



    /* Local variables */

    MDouble dcmax, x1, x2, x3, dc;

    x2 = x3 = 0;

    static MInt nadd, i, k, l, m, ikode, ic, nc, ntrail, lbl;

    static MString chedm;

    static MString colhdl[6], colhdu[6], cx2, cx3, cheval;



    if (fNu == 0) {

	fPrintf(" THERE ARE CURRENTLY NO PARAMETERS DEFINED");

	return;

    }

//*-*-                 get value of IKODE based in INKODE, ISW(2)

    ikode = inkode;

    if (inkode == 5) {

	ikode = fISW[1] + 1;

	if (ikode > 3) ikode = 3;

    }

//*-*-                 set 'default' column headings

    for (k = 1; k <= 6; ++k) {

	colhdu[k-1] = "UNDEFINED";

	colhdl[k-1] = "COLUMN HEAD";

    }

//*-*-             print title if Minos errors, and title exists.

    if (ikode == 4 && fCtitl != fCundef) {

	fPrintf(" MINUIT TASK: %s",(const char*)fCtitl);

    }

//*-*-             report function value and status

    if (fval == fUndefi) cheval = " unknown       ";

    else                 cheval = Form("%g",fval);



    if (fEDM == fBigedm) chedm = " unknown  ";

    else                 chedm = Form("%g",fEDM);



    nc = fNfcn - fNfcnfr;

    fPrintf(" FCN=%s FROM %8s  STATUS=%10s %6d CALLS   %9d TOTAL"

               ,(const char*)cheval

               ,(const char*)fCfrom

               ,(const char*)fCstatu,nc,fNfcn);

    m = fISW[1];

    if (m == 0 || m == 2 || fDcovar == 0) {

	fPrintf("                     EDM=%s    STRATEGY=%2d      %s"

                      ,(const char*)chedm,fIstrat

                      ,(const char*)fCovmes[m]);

    } else {

	dcmax = 1;

	dc    = TMath_Min(fDcovar,dcmax)*100;

	fPrintf("                     EDM=%s    STRATEGY=%2d  ERROR MATRIX UNCERTAINTY %5.1f per cent"

                      ,(const char*)chedm,fIstrat,dc);

    }



    if (ikode == 0) return;

//*-*-              find longest name (for Rene!)

    ntrail = 10;

    for (i = 1; i <= fNu; ++i) {

	if (fNvarl[i-1] < 0) continue;

	for (ic = 10; ic >= 1; --ic) {

	    if (fCpnam[i-1](ic-1,1) != " ") goto L16;

	}

	ic = 1;

L16:

	lbl = 10 - ic;

	if (lbl < ntrail) ntrail = lbl;

    }

    nadd = ntrail / 2 + 1;

    if (ikode == 1) {

	colhdu[0] = "              ";

	colhdl[0] = "      ERROR   ";

	colhdu[1] = "      PHYSICAL";

	colhdu[2] = " LIMITS       ";

	colhdl[1] = "    NEGATIVE  ";

	colhdl[2] = "    POSITIVE  ";

    }

    if (ikode == 2) {

	colhdu[0] = "              ";

	colhdl[0] = "      ERROR   ";

	colhdu[1] = "    INTERNAL  ";

	colhdl[1] = "    STEP SIZE ";

	colhdu[2] = "    INTERNAL  ";

	colhdl[2] = "      VALUE   ";

    }

    if (ikode == 3) {

	colhdu[0] = "              ";

	colhdl[0] = "      ERROR   ";

	colhdu[1] = "       STEP   ";

	colhdl[1] = "       SIZE   ";

	colhdu[2] = "      FIRST   ";

	colhdl[2] = "   DERIVATIVE ";

    }

    if (ikode == 4) {

	colhdu[0] = "    PARABOLIC ";

	colhdl[0] = "      ERROR   ";

	colhdu[1] = "        MINOS ";

	colhdu[2] = "ERRORS        ";

	colhdl[1] = "   NEGATIVE   ";

	colhdl[2] = "   POSITIVE   ";

    }



    if (ikode != 4) {

	if (fISW[1] < 3) colhdu[0] = "  APPROXIMATE ";

	if (fISW[1] < 1) colhdu[0] = " CURRENT GUESS";

    }

    fPrintf("  EXT PARAMETER              %-14s%-14s%-14s",(const char*)colhdu[0]

                                                    ,(const char*)colhdu[1]

                                                    ,(const char*)colhdu[2]);

    fPrintf("  NO.   NAME      VALUE      %-14s%-14s%-14s",(const char*)colhdl[0]

                                                    ,(const char*)colhdl[1]

                                                    ,(const char*)colhdl[2]);

//*-*-                                       . . . loop over parameters . .

    for (i = 1; i <= fNu; ++i) {

	if (fNvarl[i-1] < 0) continue;

	l = fNiofex[i-1];

        cnambf = cblank(0,nadd);

        cnambf += fCpnam[i-1];

	if (l == 0) goto L55;

//*-*-             variable parameter.

	x1  = fWerr[l-1];

	cx2 = "PLEASE GET X..";

	cx3 = "PLEASE GET X..";

	if (ikode == 1) {

	    if (fNvarl[i-1] <= 1) {

                fPrintf("%4d %-11s%14.5e%14.5e",i,(const char*)cnambf,fU[i-1],x1);

		continue;

	    } else {

		x2 = fAlim[i-1];

		x3 = fBlim[i-1];

	    }

	}

	if (ikode == 2) {

	    x2 = fDirin[l-1];

	    x3 = fX[l-1];

	}

	if (ikode == 3) {

	    x2 = fDirin[l-1];

	    x3 = fGrd[l-1];

	    if (fNvarl[i-1] > 1 && TMath_Abs(TMath_Cos(fX[l-1])) < .001) {

		cx3 = "** at limit **";

	    }

	}

	if (ikode == 4) {

	    x2 = fErn[l-1];

	    if (x2 == 0)        cx2 = " ";

	    if (x2 == fUndefi)  cx2 = "   at limit   ";

	    x3 = fErp[l-1];

	    if (x3 == 0)        cx3 = " ";

	    if (x3 == fUndefi) 	cx3 = "   at limit   ";

	}

	if (cx2 == "PLEASE GET X..")  cx2 = Form("%14.5e",x2);

	if (cx3 == "PLEASE GET X..")  cx3 = Form("%14.5e",x3);

        fPrintf("%4d %-11s%14.5e%14.5e%-14s%-14s",i

                   ,(const char*)cnambf,fU[i-1],x1

                   ,(const char*)cx2,(const char*)cx3);



//*-*-              check if parameter is at limit

	if (fNvarl[i-1] <= 1 || ikode == 3) continue;

	if (TMath_Abs(TMath_Cos(fX[l-1])) < .001) {

	    fPrintf("                                 WARNING -   - ABOVE PARAMETER IS AT LIMIT.");

	}

	continue;



//*-*-                               print constant or fixed parameter.

L55:

	colhdu[0] = "   constant   ";

	if (fNvarl[i-1] > 0)  colhdu[0] = "     fixed    ";

	if (fNvarl[i-1] == 4 && ikode == 1) {

	    fPrintf("%4d %-11s%14.5e%-14s%14.5e%14.5e",i

              ,(const char*)cnambf,fU[i-1]

              ,(const char*)colhdu[0],fAlim[i-1],fBlim[i-1]);

	} else {

	    fPrintf("%4d %-11s%14.5e%s",i

                   ,(const char*)cnambf,fU[i-1],(const char*)colhdu[0]);

	}

    }



    if (fUp != fUpdflt) {

        fPrintf("                               ERR DEF= %g",fUp);

    }

    return;

} /* mnprin_ */



//______________________________________________________________________________

void Midnight::mnpsdf()

{

//*-*-*-*-*-*Calculates the eigenvalues of v to see if positive-def*-*-*-*-*

//*-*        ======================================================

//*-*        if not, adds constant along diagonal to make positive.

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    /* Local variables */

    static MDouble s[kMAXDIM], dgmin, padd, pmin, pmax, dg, epspdf, epsmin;

    static MInt ndex, i, j, ndexd, ip, ifault;

    static MString chbuff, ctemp;



    epsmin = 1e-6;

    epspdf = TMath_Max(epsmin,fEpsma2);

    dgmin  = fVhmat[0];

//*-*-                       Check if negative or zero on diagonal

    for (i = 1; i <= fNpar; ++i) {

	ndex = i*(i + 1) / 2;

	if (fVhmat[ndex-1] <= 0) {

	    mnwarn("W", fCfrom, Form("Negative diagonal element %d in Error Matrix",i));

	}

	if (fVhmat[ndex-1] < dgmin) dgmin = fVhmat[ndex-1];

    }

    if (dgmin <= 0) {

	dg    = epspdf + 1 - dgmin;

	mnwarn("W", fCfrom, Form("%g added to diagonal of error matrix",dg));

    } else {

	dg = 0;

    }

//*-*-                   Store VHMAT in P, make sure diagonal pos.

    for (i = 1; i <= fNpar; ++i) {

	ndex  = i*(i-1) / 2;

	ndexd = ndex + i;

	fVhmat[ndexd-1] += dg;

	s[i-1] = 1 / TMath_Sqrt(fVhmat[ndexd-1]);

	for (j = 1; j <= i; ++j) {

	    ++ndex;

	    fP[i + j*fMaxpar - fMaxpar-1] = fVhmat[ndex-1]*s[i-1]*s[j-1];

	}

    }

//*-*-     call eigen (p,p,maxint,npar,pstar,-npar)

    mneig(fP, fMaxint, fNpar, fMaxint, fPstar, epspdf, ifault);

    pmin = fPstar[0];

    pmax = fPstar[0];

    for (ip = 2; ip <= fNpar; ++ip) {

	if (fPstar[ip-1] < pmin) pmin = fPstar[ip-1];

	if (fPstar[ip-1] > pmax) pmax = fPstar[ip-1];

    }

    pmax = TMath_Max(TMath_Abs(pmax),MDouble(1));

    if (pmin <= 0 && fLwarn || fISW[4] >= 2) {

	fPrintf(" EIGENVALUES OF SECOND-DERIVATIVE MATRIX:");

        ctemp = "       ";

	for (ip = 1; ip <= fNpar; ++ip) {

	    ctemp += Form(" %11.4e",fPstar[ip-1]);

	}

        fPrintf((const char*)ctemp);

    }

    if (pmin > epspdf*pmax) return;

    if (fISW[1] == 3) fISW[1] = 2;

    padd = pmax*.001 - pmin;

    for (ip = 1; ip <= fNpar; ++ip) {

	ndex = ip*(ip + 1) / 2;

	fVhmat[ndex-1] *= padd + 1;

    }

    fCstatu = "NOT POSDEF";

    mnwarn("W", fCfrom, Form("MATRIX FORCED POS-DEF BY ADDING %f TO DIAGONAL.",padd));



} /* mnpsdf_ */



//______________________________________________________________________________

void Midnight::mnrazz(MDouble ynew, MDouble *pnew, MDouble *y, MInt &jh, MInt &jl)

{

//*-*-*-*-*Called only by MNSIMP (and MNIMPR) to add a new point*-*-*-*-*-*-*

//*-*      =====================================================

//*-*        and remove an old one from the current simplex, and get the

//*-*        estimated distance to minimum.

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    /* Local variables */

    static MDouble pbig, plit;

    static MInt i, j, nparp1;



    /* Function Body */

    for (i = 1; i <= fNpar; ++i) { fP[i + jh*fMaxpar - fMaxpar-1] = pnew[i-1]; }

    y[jh-1] = ynew;

    if (ynew < fAmin) {

	for (i = 1; i <= fNpar; ++i) { fX[i-1] = pnew[i-1]; }

	mninex(fX);

	fAmin   = ynew;

	fCstatu = "PROGRESS  ";

	jl      = jh;

    }

    jh     = 1;

    nparp1 = fNpar + 1;

    for (j = 2; j <= nparp1; ++j) { if (y[j-1] > y[jh-1]) jh = j; }

    fEDM = y[jh-1] - y[jl-1];

    if (fEDM <= 0) goto L45;

    for (i = 1; i <= fNpar; ++i) {

	pbig = fP[i-1];

	plit = pbig;

	for (j = 2; j <= nparp1; ++j) {

	    if (fP[i + j*fMaxpar - fMaxpar-1] > pbig) pbig = fP[i + j*fMaxpar - fMaxpar-1];

	    if (fP[i + j*fMaxpar - fMaxpar-1] < plit) plit = fP[i + j*fMaxpar - fMaxpar-1];

	}

	fDirin[i-1] = pbig - plit;

    }

L40:

    return;

L45:

    fPrintf("   FUNCTION VALUE DOES NOT SEEM TO DEPEND ON ANY OF THE%d VARIABLE PARAMETERS.",fNpar);

    fPrintf("          VERIFY THAT STEP SIZES ARE BIG ENOUGH AND CHECK FCN LOGIC.");

    fPrintf(" *******************************************************************************");

    fPrintf(" *******************************************************************************");

    goto L40;

} /* mnrazz_ */



//______________________________________________________________________________

void Midnight::mnrn15(MDouble &val, MInt &inseed)

{

//*-*-*-*-*-*-*This is a super-portable random number generator*-*-*-*-*-*-*

//*-*          ================================================

//*-*         It should not overflow on any 32-bit machine.

//*-*         The cycle is only ~10**9, so use with care!

//*-*         Note especially that VAL must not be undefined on input.

//*-*                    Set Default Starting Seed

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    /* Initialized data */



    static MInt iseed = 12345;



    MInt k;



    if (val == 3) goto L100;

    inseed = iseed;

    k      = iseed / 53668;

    iseed  = (iseed - k*53668)*40014 - k*12211;

    if (iseed < 0) iseed += 2147483563;

    val = MDouble(iseed*4.656613e-10);

    return;

//*-*               "entry" to set seed, flag is VAL=3

L100:

    iseed = inseed;

} /* mnrn15_ */



//______________________________________________________________________________

void Midnight::mnrset(MInt iopt)

{

//*-*-*-*-*-*-*-*Resets function value and errors to UNDEFINED*-*-*-*-*-*-*-*

//*-*            =============================================

//*-*    If IOPT=1,

//*-*    If IOPT=0, sets only MINOS errors to undefined

//*-*        Called from MNCLER and whenever problem changes, for example

//*-*        after SET LIMITS, SET PARAM, CALL FCN 6

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    static MInt iext, i;



    fCstatu = "RESET     ";

    if (iopt >= 1) {

	fAmin   = fUndefi;

	fFval3  = TMath_Abs(fAmin)*2 + 1;

	fEDM    = fBigedm;

	fISW[3] = 0;

	fISW[1] = 0;

	fDcovar = 1;

	fISW[0] = 0;

    }

    fLnolim = kTRUE;

    for (i = 1; i <= fNpar; ++i) {

	iext = fNexofi[i-1];

	if (fNvarl[iext-1] >= 4) fLnolim = kFALSE;

	fErp[i-1] = 0;

	fErn[i-1] = 0;

	fGlobcc[i-1] = 0;

    }

    if (fISW[1] >= 1) {

	fISW[1] = 1;

	fDcovar = TMath_Max(fDcovar,.5);

    }

} /* mnrset_ */



//______________________________________________________________________________

void Midnight::mnsave()

{

//*-*-*-*Writes current parameter values and step sizes onto file ISYSSA*-*-*

//*-*    ===============================================================

//*-*          in format which can be reread by Minuit for restarting.

//*-*       The covariance matrix is also output if it exists.

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



   fPrintf("mnsave is dummy in the base class Minuit: Use MinuitOld");



} /* mnsave_ */



//______________________________________________________________________________

void Midnight::mnscan()

{

//*-*-*-*-*Scans the values of FCN as a function of one parameter*-*-*-*-*-*

//*-*      ======================================================

//*-*        and plots the resulting values as a curve using MNPLOT.

//*-*        It may be called to scan one parameter or all parameters.

//*-*        retains the best function and parameter values found.

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    /* Local variables */

    static MDouble step, uhigh, xhreq, xlreq, ubest, fnext, unext, xh, xl;

    static MInt ipar, iint, icall, ncall, nbins, nparx;

    static MInt nxypt, nccall, iparwd;



    xlreq = TMath_Min(fWord7[2],fWord7[3]);

    xhreq = TMath_Max(fWord7[2],fWord7[3]);

    ncall = MInt((fWord7[1] + .01));

    if (ncall <= 1)  ncall = 41;

    if (ncall > 101) ncall = 101;

    nccall = ncall;

    if (fAmin == fUndefi) mnamin();

    iparwd  = MInt((fWord7[0] + .1));

    ipar    = TMath_Max(iparwd,0);

    iint    = fNiofex[ipar-1];

    fCstatu = "NO CHANGE";

    if (iparwd > 0) goto L200;



//*-*-        equivalent to a loop over parameters requested

L100:

    ++ipar;

    if (ipar > fNu) goto L900;

    iint = fNiofex[ipar-1];

    if (iint <= 0) goto L100;

//*-*-        set up range for parameter IPAR

L200:

    ubest    = fU[ipar-1];

    fXpt[0]  = ubest;

    fYpt[0]  = fAmin;

    fChpt[0] = ' ';

    fXpt[1]  = ubest;

    fYpt[1]  = fAmin;

    fChpt[1] = 'X';

    nxypt    = 2;

    if (fNvarl[ipar-1] > 1) goto L300;



//*-*-        no limits on parameter

    if (xlreq == xhreq) goto L250;

    unext = xlreq;

    step = (xhreq - xlreq) / MDouble(ncall-1);

    goto L500;

L250:

    xl = ubest - fWerr[iint-1];

    xh = ubest + fWerr[iint-1];

    mnbins(xl, xh, ncall, unext, uhigh, nbins, step);

    nccall = nbins + 1;

    goto L500;

//*-*-        limits on parameter

L300:

    if (xlreq == xhreq) goto L350;

//*-*  Computing MAX

    xl = TMath_Max(xlreq,fAlim[ipar-1]);

//*-*  Computing MIN

    xh = TMath_Min(xhreq,fBlim[ipar-1]);

    if (xl >= xh) goto L700;

    unext = xl;

    step  = (xh - xl) / MDouble(ncall-1);

    goto L500;

L350:

    unext = fAlim[ipar-1];

    step = (fBlim[ipar-1] - fAlim[ipar-1]) / MDouble(ncall-1);

//*-*-        main scanning loop over parameter IPAR

L500:

    for (icall = 1; icall <= nccall; ++icall) {

	fU[ipar-1] = unext;

	nparx = fNpar;

	(*fFCN)(nparx, fGin, fnext, fU, 4);	++fNfcn;

	++nxypt;

	fXpt[nxypt-1]  = unext;

	fYpt[nxypt-1]  = fnext;

	fChpt[nxypt-1] = '*';

	if (fnext < fAmin) {

	    fAmin   = fnext;

	    ubest   = unext;

	    fCstatu = "IMPROVED  ";

	}

	unext += step;

    }

//*-*-        finished with scan of parameter IPAR

    fU[ipar-1] = ubest;

    mnexin(fX);

    fPrintf("%dSCAN OF PARAMETER NO. %d,  %s"

         ,fNewpag,ipar,(const char*)fCpnam[ipar-1]);

    mnplot(fXpt, fYpt, fChpt, nxypt, fNpagwd, fNpagln);

    goto L800;

L700:

    fPrintf(" REQUESTED RANGE OUTSIDE LIMITS FOR PARAMETER  %d",ipar);

L800:

    if (iparwd <= 0) goto L100;

//*-*-        finished with all parameters

L900:

    mnprin(5, fAmin);

} /* mnscan_ */



//______________________________________________________________________________

void Midnight::mnseek()

{

//*-*-*-*Performs a rough (but global) minimization by monte carlo search*-*

//*-*    ================================================================

//*-*        Each time a new minimum is found, the search area is shifted

//*-*        to be centered at the best value.  Random points are chosen

//*-*        uniformly over a hypercube determined by current step sizes.

//*-*   The Metropolis algorithm accepts a worse point with probability

//*-*      exp(-d/UP), where d is the degradation.  Improved points

//*-*      are of course always accepted.  Actual steps are random

//*-*      multiples of the nominal steps (DIRIN).

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    /* Local variables */

    static MDouble xmid[kMAXDIM], dxdi, rnum, ftry, rnum1, rnum2, alpha;

    static MDouble flast, xbest[kMAXDIM], bar;

    static MInt ipar, iext, j, ifail, iseed, nparx, istep, ib, mxfail, mxstep;



    mxfail = MInt(fWord7[0]);

    if (mxfail <= 0) mxfail = fNpar*20 + 100;

    mxstep = mxfail*10;

    if (fAmin == fUndefi) mnamin();

    alpha = fWord7[1];

    if (alpha <= 0) alpha = 3;

    if (fISW[4] >= 1) {

	fPrintf(" MNSEEK: MONTE CARLO MINIMIZATION USING METROPOLIS ALGORITHM");

        fPrintf(" TO STOP AFTER %6d SUCCESSIVE FAILURES, OR %7d STEPS",mxfail,mxstep);

        fPrintf(" MAXIMUM STEP SIZE IS %9.3f ERROR BARS.",alpha);

    }

    fCstatu = "INITIAL  ";

    if (fISW[4] >= 2) mnprin(2, fAmin);

    fCstatu = "UNCHANGED ";

    ifail   = 0;

    rnum    = 0;

    rnum1   = 0;

    rnum2   = 0;

    nparx   = fNpar;

    flast   = fAmin;

//*-*-             set up step sizes, starting values

    for (ipar = 1; ipar <= fNpar; ++ipar) {

	iext = fNexofi[ipar-1];

	fDirin[ipar-1] = alpha*2*fWerr[ipar-1];

	if (fNvarl[iext-1] > 1) {

//*-*-             parameter with limits

	    mndxdi(fX[ipar-1], ipar-1, dxdi);

	    if (dxdi == 0) dxdi = 1;

	    fDirin[ipar-1] = alpha*2*fWerr[ipar-1] / dxdi;

	    if (TMath_Abs(fDirin[ipar-1]) > 6.2831859999999997) {

		fDirin[ipar-1] = 6.2831859999999997;

	    }

	}

	xmid[ipar-1] = fX[ipar-1];

	xbest[ipar-1] = fX[ipar-1];

    }

//*-*-                             search loop

    for (istep = 1; istep <= mxstep; ++istep) {

	if (ifail >= mxfail) break;

	for (ipar = 1; ipar <= fNpar; ++ipar) {

	    mnrn15(rnum1, iseed);

	    mnrn15(rnum2, iseed);

	    fX[ipar-1] = xmid[ipar-1] + (rnum1 + rnum2 - 1)*.5*fDirin[ipar-1];

	}

	mninex(fX);

	(*fFCN)(nparx, fGin, ftry, fU, 4);	++fNfcn;

	if (ftry < flast) {

	    if (ftry < fAmin) {

		fCstatu = "IMPROVEMNT";

		fAmin = ftry;

		for (ib = 1; ib <= fNpar; ++ib) { xbest[ib-1] = fX[ib-1]; }

		ifail = 0;

		if (fISW[4] >= 2) mnprin(2, fAmin);

	    }

	    goto L300;

	} else {

	    ++ifail;

//*-*-                  Metropolis algorithm

	    bar = (fAmin - ftry) / fUp;

	    mnrn15(rnum, iseed);

	    if (bar < TMath_Log(rnum)) continue;

	}

//*-*-                   Accept new point, move there

L300:

	for (j = 1; j <= fNpar; ++j) { xmid[j-1] = fX[j-1];	}

	flast = ftry;

    }

//*-*-                              end search loop

    if (fISW[4] > 1) {

        fPrintf(" MNSEEK: %5d SUCCESSIVE UNSUCCESSFUL TRIALS.",ifail);

    }

    for (ib = 1; ib <= fNpar; ++ib) { fX[ib-1] = xbest[ib-1]; }

    mninex(fX);

    if (fISW[4] >= 1) mnprin(2, fAmin);

    if (fISW[4] == 0) mnprin(0, fAmin);

} /* mnseek_ */



//______________________________________________________________________________

void Midnight::mnset()

{

//*-*-*-*-*Interprets the commands that start with SET and SHOW*-*-*-*-*-*-*

//*-*      ====================================================

//*-*        Called from MNEXCM

//*-*        file characteristics for SET INPUT

//*-*       'SET ' or 'SHOW',  'ON ' or 'OFF', 'SUPPRESSED' or 'REPORTED  '

//*-*        explanation of print level numbers -1:3  and strategies 0:2

//*-*        identification of debug options

//*-*        things that can be set or shown

//*-*        options not intended for normal users

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    /* Initialized data */



    static MString cname[30] = {

      "FCN value ",

      "PARameters",

      "LIMits    ",

      "COVariance",

      "CORrelatio",

      "PRInt levl",

      "NOGradient",

      "GRAdient  ",

      "ERRor def ",

      "INPut file",

      "WIDth page",

      "LINes page",

      "NOWarnings",

      "WARnings  ",

      "RANdom gen",

      "TITle     ",

      "STRategy  ",

      "EIGenvalue",

      "PAGe throw",

      "MINos errs",

      "EPSmachine",

      "OUTputfile",

      "BATch     ",

      "INTeractiv",

      "VERsion   ",

      "reserve   ",

      "NODebug   ",

      "DEBug     ",

      "SHOw      ",

      "SET       "};



    static MInt nname = 25;

    static MInt nntot = 30;

    static MString cprlev[5] = {

      "-1: NO OUTPUT EXCEPT FROM SHOW    ",

      " 0: REDUCED OUTPUT                ",

      " 1: NORMAL OUTPUT                 ",

      " 2: EXTRA OUTPUT FOR PROBLEM CASES",

      " 3: MAXIMUM OUTPUT                "};



    static MString cstrat[3] = {

      " 0: MINIMIZE THE NUMBER OF CALLS TO FUNCTION",

      " 1: TRY TO BALANCE SPEED AGAINST RELIABILITY",

      " 2: MAKE SURE MINIMUM TRUE, ERRORS CORRECT  "};



    static MString cdbopt[7] = {

      "REPORT ALL EXCEPTIONAL CONDITIONS      ",

      "MNLINE: LINE SEARCH MINIMIZATION       ",

      "MNDERI: FIRST DERIVATIVE CALCULATIONS  ",

      "MNHESS: SECOND DERIVATIVE CALCULATIONS ",

      "MNMIGR: COVARIANCE MATRIX UPDATES      ",

      "MNHES1: FIRST DERIVATIVE UNCERTAINTIES ",

      "MNCONT: MNCONTOUR PLOT (MNCROS SEARCH) "};



    /* System generated locals */

    MInt f_inqu();



    /* Local variables */

    static MDouble val;

    static MInt iset, iprm, i, jseed, kname, iseed, iunit, id, ii, kk;

    static MInt ikseed, idbopt, igrain, iswsav, isw2;

    static MString  cfname, cmode, ckind,  cwarn, copt, ctemp, ctemp2;

    static MBool lname;



    for (i = 1; i <= nntot; ++i) {

	ctemp  = cname[i-1](0,3);

	ctemp2 = fCword(3,7);

	if (strstr((const char*)ctemp2, (const char*)ctemp)) goto L5;

    }

    i = 0;

L5:

    kname = i;



//*-*-          Command could be SET xxx, SHOW xxx,  HELP SET or HELP SHOW

    if (fCword(0,3) == "HEL")  goto L2000;

    if (fCword(0,3) == "SHO")  goto L1000;

    if (fCword(0,3) != "SET")  goto L1900;

//*-*-                          ---

    ckind = "SET ";

//*-*-                                       . . . . . . . . . . set unknown

    if (kname <= 0) goto L1900;

//*-*-                                       . . . . . . . . . . set known

    switch ((int)kname) {

	case 1:  goto L3000;

	case 2:  goto L20;

	case 3:  goto L30;

	case 4:  goto L40;

	case 5:  goto L3000;

	case 6:  goto L60;

	case 7:  goto L70;

	case 8:  goto L80;

	case 9:  goto L90;

	case 10:  goto L100;

	case 11:  goto L110;

	case 12:  goto L120;

	case 13:  goto L130;

	case 14:  goto L140;

	case 15:  goto L150;

	case 16:  goto L160;

	case 17:  goto L170;

	case 18:  goto L3000;

	case 19:  goto L190;

	case 20:  goto L3000;

	case 21:  goto L210;

	case 22:  goto L220;

	case 23:  goto L230;

	case 24:  goto L240;

	case 25:  goto L3000;

	case 26:  goto L1900;

	case 27:  goto L270;

	case 28:  goto L280;

	case 29:  goto L290;

	case 30:  goto L300;

    }



//*-*-                                       . . . . . . . . . . set param

L20:

    iprm = MInt(fWord7[0]);

    if (iprm > fNu) goto L25;

    if (iprm <= 0) goto L25;

    if (fNvarl[iprm-1] < 0) goto L25;

    fU[iprm-1] = fWord7[1];

    mnexin(fX);

    isw2 = fISW[1];

    mnrset(1);

//*-*-       Keep approximate covariance matrix, even if new param value

    fISW[1] = TMath_Min(isw2,1);

    fCfrom  = "SET PARM";

    fNfcnfr = fNfcn;

    fCstatu = "NEW VALUES";

    return;

L25:

    fPrintf(" UNDEFINED PARAMETER NUMBER.  IGNORED.");

    return;

//*-*-                                       . . . . . . . . . . set limits

L30:

    mnlims();

    return;

//*-*-                                       . . . . . . . . . . set covar

L40:

//*-*   this command must be handled by MNREAD, and is not Fortran-callable

    goto L3000;

//*-*-                                       . . . . . . . . . . set print

L60:

    fISW[4] = MInt(fWord7[0]);

    return;

//*-*-                                       . . . . . . . . . . set nograd

L70:

    fISW[2] = 0;

    return;

//*-*-                                       . . . . . . . . . . set grad

L80:

    mngrad();

    return;

//*-*-                                       . . . . . . . . . . set errdef

L90:

    if (fWord7[0] == fUp) return;

    if (fWord7[0] <= 0) {

	if (fUp == fUpdflt) return;

	fUp = fUpdflt;

    } else {

	fUp = fWord7[0];

    }

    for (i = 1; i <= fNpar; ++i) {

	fErn[i-1] = 0;

	fErp[i-1] = 0;

    }

    mnwerr();

    return;

//*-*-                                       . . . . . . . . . . set input

//*-* This command must be handled by MNREAD. If it gets this far,

//*-*-        it is illegal.

L100:

    goto L3000;

//*-*-                                       . . . . . . . . . . set width

L110:

    fNpagwd = MInt(fWord7[0]);

    fNpagwd = TMath_Max(fNpagwd,50);

    return;



L120:

    fNpagln = MInt(fWord7[0]);

    return;

//*-*-                                       . . . . . . . . . . set nowarn



L130:

    fLwarn = kFALSE;

    return;

//*-*-                                       . . . . . . . . . . set warn

L140:

    fLwarn = kTRUE;

    mnwarn("W", "SHO", "SHO");

    return;

//*-*-                                       . . . . . . . . . . set random

L150:

    jseed = MInt(fWord7[0]);

    val = 3;

    mnrn15(val, jseed);

    if (fISW[4] > 0) {

        fPrintf(" MINUIT RANDOM NUMBER SEED SET TO %d",jseed);

    }

    return;

//*-*-                                       . . . . . . . . . . set title

L160:

//*-*   this command must be handled by MNREAD, and is not Fortran-callable

    goto L3000;

//*-*-                                       . . . . . . . . . set strategy

L170:

    fIstrat = MInt(fWord7[0]);

    fIstrat = TMath_Max(fIstrat,0);

    fIstrat = TMath_Min(fIstrat,2);

    if (fISW[4] > 0) goto L1172;

    return;

//*-*-                                      . . . . . . . . . set page throw

L190:

    fNewpag = MInt(fWord7[0]);

    goto L1190;

//*-*-                                       . . . . . . . . . . set epsmac

L210:

    if (fWord7[0] > 0 && fWord7[0] < .1) {

	fEpsmac = fWord7[0];

    }

    fEpsma2 = TMath_Sqrt(fEpsmac);

    goto L1210;

//*-*-                                      . . . . . . . . . . set outputfile

L220:

    iunit = MInt(fWord7[0]);

    fIsyswr = iunit;

    fIstkwr[0] = iunit;

    if (fISW[4] >= 0) goto L1220;

    return;

//*-*-                                       . . . . . . . . . . set batch

L230:

    fISW[5] = 0;

    if (fISW[4] >= 0) goto L1100;

    return;

//*-*-                                      . . . . . . . . . . set interactive

L240:

    fISW[5] = 1;

    if (fISW[4] >= 0) goto L1100;

    return;

//*-*-                                       . . . . . . . . . . set nodebug

L270:

    iset = 0;

    goto L281;

//*-*-                                       . . . . . . . . . . set debug

L280:

    iset = 1;

L281:

    idbopt = MInt(fWord7[0]);

    if (idbopt > 6) goto L288;

    if (idbopt >= 0) {

	fIdbg[idbopt] = iset;

	if (iset == 1) fIdbg[0] = 1;

    } else {

//*-*-            SET DEBUG -1  sets all debug options

	for (id = 0; id <= 6; ++id) { fIdbg[id] = iset; }

    }

    fLrepor = fIdbg[0] >= 1;

    mnwarn("D", "SHO", "SHO");

    return;

L288:

    fPrintf(" UNKNOWN DEBUG OPTION %d REQUESTED. IGNORED",idbopt);

    return;

//*-*-                                       . . . . . . . . . . set show

L290:

//*-*-                                       . . . . . . . . . . set set

L300:

    goto L3000;

//*-*-               -----------------------------------------------------

L1000:

//*-*-              at this point, CWORD must be 'SHOW'

    ckind = "SHOW";

    if (kname <= 0) goto L1900;



    switch ((int)kname) {

	case 1:  goto L1010;

	case 2:  goto L1020;

	case 3:  goto L1030;

	case 4:  goto L1040;

	case 5:  goto L1050;

	case 6:  goto L1060;

	case 7:  goto L1070;

	case 8:  goto L1070;

	case 9:  goto L1090;

	case 10:  goto L1100;

	case 11:  goto L1110;

	case 12:  goto L1120;

	case 13:  goto L1130;

	case 14:  goto L1130;

	case 15:  goto L1150;

	case 16:  goto L1160;

	case 17:  goto L1170;

	case 18:  goto L1180;

	case 19:  goto L1190;

	case 20:  goto L1200;

	case 21:  goto L1210;

	case 22:  goto L1220;

	case 23:  goto L1100;

	case 24:  goto L1100;

	case 25:  goto L1250;

	case 26:  goto L1900;

	case 27:  goto L1270;

	case 28:  goto L1270;

	case 29:  goto L1290;

	case 30:  goto L1300;

    }



//*-*-                                       . . . . . . . . . . show fcn

L1010:

    if (fAmin == fUndefi) mnamin();

    mnprin(0, fAmin);

    return;

//*-*-                                       . . . . . . . . . . show param

L1020:

    if (fAmin == fUndefi) mnamin();

    mnprin(5, fAmin);

    return;

//*-*-                                       . . . . . . . . . . show limits

L1030:

    if (fAmin == fUndefi) mnamin();

    mnprin(1, fAmin);

    return;

//*-*-                                       . . . . . . . . . . show covar

L1040:

    mnmatu(1);

    return;

//*-*-                                       . . . . . . . . . . show corre

L1050:

    mnmatu(0);

    return;

//*-*-                                       . . . . . . . . . . show print

L1060:

    if (fISW[4] < -1) fISW[4] = -1;

    if (fISW[4] > 3)  fISW[4] = 3;

    fPrintf(" ALLOWED PRINT LEVELS ARE:");

    fPrintf("                           %s",(const char*)cprlev[0]);

    fPrintf("                           %s",(const char*)cprlev[1]);

    fPrintf("                           %s",(const char*)cprlev[2]);

    fPrintf("                           %s",(const char*)cprlev[3]);

    fPrintf("                           %s",(const char*)cprlev[4]);

    fPrintf(" CURRENT PRINTOUT LEVEL IS %s",(const char*)cprlev[fISW[4]]);

    return;

//*-*-                                       . . . . . . . show nograd, grad

L1070:

    if (fISW[2] <= 0) {

	fPrintf(" NOGRAD IS SET.  DERIVATIVES NOT COMPUTED IN FCN.");

    } else {

	fPrintf("   GRAD IS SET.  USER COMPUTES DERIVATIVES IN FCN.");

    }

    return;

//*-*-                                      . . . . . . . . . . show errdef

L1090:

    fPrintf(" ERRORS CORRESPOND TO FUNCTION CHANGE OF %g",fUp);

    return;

//*-*-                                      . . . . . . . . . . show input,

//*-*-                                               batch, or interactive

L1100:

//    ioin__1.inerr = 0;

//    ioin__1.inunit = fIsysrd;

//    ioin__1.infile = 0;

//    ioin__1.inex = 0;

//    ioin__1.inopen = 0;

//    ioin__1.innum = 0;

//    ioin__1.innamed = &lname;

//    ioin__1.innamlen = 64;

//    ioin__1.inname = cfname;

//    ioin__1.inacc = 0;

//    ioin__1.inseq = 0;

//    ioin__1.indir = 0;

//    ioin__1.infmt = 0;

//    ioin__1.inform = 0;

//    ioin__1.inunf = 0;

//    ioin__1.inrecl = 0;

//    ioin__1.innrec = 0;

//    ioin__1.inblank = 0;

//    f_inqu(&ioin__1);

    cmode = "BATCH MODE      ";

    if (fISW[5] == 1) cmode  = "INTERACTIVE MODE";

    if (! lname)      cfname = "unknown";

    fPrintf(" INPUT NOW BEING READ IN %s FROM UNIT NO. %d FILENAME: %s"

           ,(const char*)cmode,fIsysrd,(const char*)cfname);

    return;

//*-*-                                      . . . . . . . . . . show width

L1110:

    fPrintf("          PAGE WIDTH IS SET TO %d COLUMNS",fNpagwd);

    return;

//*-*-                                      . . . . . . . . . . show lines

L1120:

    fPrintf("          PAGE LENGTH IS SET TO %d LINES",fNpagln);

    return;

//*-*-                                      . . . . . . .show nowarn, warn

L1130:

    cwarn = "SUPPRESSED";

    if (fLwarn) cwarn = "REPORTED  ";

    fPrintf("%s",(const char*)cwarn);

    if (! fLwarn) mnwarn("W", "SHO", "SHO");

    return;

//*-*-                                     . . . . . . . . . . show random

L1150:

    val = 0;

    mnrn15(val, igrain);

    ikseed = igrain;

    fPrintf(" MINUIT RNDM SEED IS CURRENTLY=",ikseed);

    val   = 3;

    iseed = ikseed;

    mnrn15(val, iseed);

    return;

//*-*-                                       . . . . . . . . . show title

L1160:

    fPrintf(" TITLE OF CURRENT TASK IS:%s",(const char*)fCtitl);

    return;

//*-*-                                       . . . . . . . show strategy

L1170:

    fPrintf(" ALLOWED STRATEGIES ARE:");

    fPrintf("                    %s",(const char*)cstrat[0]);

    fPrintf("                    %s",(const char*)cstrat[1]);

    fPrintf("                    %s",(const char*)cstrat[2]);

L1172:

    fPrintf(" NOW USING STRATEGY %s",(const char*)cstrat[fIstrat]);

    return;

//*-*-                                         . . . . . show eigenvalues

L1180:

    iswsav = fISW[4];

    fISW[4] = 3;

    if (fISW[1] < 1) {

	fPrintf("%s",(const char*)fCovmes[0]);

    } else {

	mnpsdf();

    }

    fISW[4] = iswsav;

    return;

//*-*-                                           . . . . . show page throw

L1190:

    fPrintf(" PAGE THROW CARRIAGE CONTROL = %d",fNewpag);

    if (fNewpag == 0) {

	fPrintf(" NO PAGE THROWS IN MINUIT OUTPUT");

    }

    return;

//*-*-                                       . . . . . . show minos errors

L1200:

    for (ii = 1; ii <= fNpar; ++ii) {

	if (fErp[ii-1] > 0 || fErn[ii-1] < 0) goto L1204;

    }

    fPrintf("       THERE ARE NO MINOS ERRORS CURRENTLY VALID.");

    return;

L1204:

    mnprin(4, fAmin);

    return;

//*-*-                                       . . . . . . . . . show epsmac

L1210:

    fPrintf(" FLOATING-POINT NUMBERS ASSUMED ACCURATE TO %g",fEpsmac);

    return;

//*-*-                                       . . . . . . show outputfiles

L1220:

    fPrintf("  MINUIT PRIMARY OUTPUT TO UNIT %d",fIsyswr);

    return;

//*-*-                                       . . . . . . show version

L1250:

    fPrintf(" THIS IS MINUIT VERSION:%s",(const char*)fCvrsn);

    return;

//*-*-                                       . . . . . . show nodebug, debug

L1270:

    for (id = 0; id <= 6; ++id) {

	copt = "OFF";

	if (fIdbg[id] >= 1) copt = "ON ";

        fPrintf("          DEBUG OPTION %3d IS %3s :%s"

             ,id,(const char*)copt,(const char*)cdbopt[id]);

    }

    if (! fLrepor) mnwarn("D", "SHO", "SHO");

    return;

//*-*-                                       . . . . . . . . . . show show

L1290:

    ckind = "SHOW";

    goto L2100;

//*-*-                                       . . . . . . . . . . show set

L1300:

    ckind = "SET ";

    goto L2100;

//*-*-               -----------------------------------------------------

//*-*-                             UNKNOWN COMMAND

L1900:

    fPrintf(" THE COMMAND:%10s IS UNKNOWN.",(const char*)fCword);

    goto L2100;

//*-*-               -----------------------------------------------------

//*-*-                   HELP SHOW,  HELP SET,  SHOW SET, or SHOW SHOW

L2000:

    ckind = "SET ";

    ctemp = fCword(3,7);

    if (strcmp((const char*)ctemp, "SHO")) ckind = "SHOW";

L2100:

    fPrintf(" THE FORMAT OF THE %4s COMMAND IS:",(const char*)ckind);

    fPrintf(" %s xxx    [numerical arguments if any]",(const char*)ckind);

    fPrintf(" WHERE xxx MAY BE ONE OF THE FOLLOWING:");

    for (kk = 1; kk <= nname; ++kk) {

        fPrintf(" %s",(const char*)cname[kk-1]);

    }

    return;

//*-*-               -----------------------------------------------------

//*-*-                              ILLEGAL COMMAND

L3000:

    fPrintf(" ABOVE COMMAND IS ILLEGAL.   IGNORED");



} /* mnset_ */



//______________________________________________________________________________

void Midnight::mnsimp()

{

//*-*-*-*-*Minimization using the simplex method of Nelder and Mead*-*-*-*-*

//*-*      ========================================================

//*-*        Performs a minimization using the simplex method of Nelder

//*-*        and Mead (ref. -- Comp. J. 7,308 (1965)).

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    /* Initialized data */



    static MDouble alpha = 1;

    static MDouble beta = .5;

    static MDouble gamma = 2;

    static MDouble rhomin = 4;

    static MDouble rhomax = 8;



    /* Local variables */

    static MDouble dmin_, dxdi, yrho, f, ynpp1, y[kMAXDIM+1], aming, ypbar;

    static MDouble bestx, ystar, y1, y2, ystst, pb, wg;

    static MDouble absmin, rho, sig2, rho1, rho2;

    static MInt npfn, i, j, k, jhold, ncycl, nparx;

    static MInt nparp1, kg, jh, nf, jl, ns;



    if (fNpar <= 0) return;

    if (fAmin == fUndefi) mnamin();

    fCfrom  = "SIMPLEX ";

    fNfcnfr = fNfcn;

    fCstatu = "UNCHANGED ";

    npfn    = fNfcn;

    nparp1  = fNpar + 1;

    nparx   = fNpar;

    rho1    = alpha + 1;

    rho2    = rho1 + alpha*gamma;

    wg      = 1 / MDouble(fNpar);

    if (fISW[4] >= 0) {

	fPrintf(" START SIMPLEX MINIMIZATION.    CONVERGENCE WHEN EDM .LT. %g",fEpsi);

    }

    for (i = 1; i <= fNpar; ++i) {

	fDirin[i-1] = fWerr[i-1];

	mndxdi(fX[i-1], i-1, dxdi);

	if (dxdi != 0) fDirin[i-1] = fWerr[i-1] / dxdi;

	dmin_ = fEpsma2*TMath_Abs(fX[i-1]);

	if (fDirin[i-1] < dmin_) fDirin[i-1] = dmin_;

    }

//*-* **       choose the initial simplex using single-parameter searches

L1:

    ynpp1 = fAmin;

    jl = nparp1;

    y[nparp1-1] = fAmin;

    absmin = fAmin;

    for (i = 1; i <= fNpar; ++i) {

	aming      = fAmin;

	fPbar[i-1] = fX[i-1];

	bestx      = fX[i-1];

	kg         = 0;

	ns         = 0;

	nf         = 0;

L4:

	fX[i-1] = bestx + fDirin[i-1];

	mninex(fX);

	(*fFCN)(nparx, fGin, f, fU, 4);	++fNfcn;

	if (f <= aming) goto L6;

//*-*-        failure

	if (kg == 1) goto L8;

	kg = -1;

	++nf;

	fDirin[i-1] *= -.4;

	if (nf < 3) goto L4;

	ns = 6;

//*-*-        success

L6:

	bestx        = fX[i-1];

	fDirin[i-1] *= 3;

	aming        = f;

	fCstatu      = "PROGRESS  ";

	kg           = 1;

	++ns;

	if (ns < 6) goto L4;

//*-*-        local minimum found in ith direction

L8:

	y[i-1] = aming;

	if (aming < absmin) jl = i;

	if (aming < absmin) absmin = aming;

	fX[i-1] = bestx;

	for (k = 1; k <= fNpar; ++k) { fP[k + i*fMaxpar - fMaxpar-1] = fX[k-1]; }

    }

    jh    = nparp1;

    fAmin = y[jl-1];

    mnrazz(ynpp1, fPbar, y, jh, jl);

    for (i = 1; i <= fNpar; ++i) { fX[i-1] = fP[i + jl*fMaxpar - fMaxpar-1]; }

    mninex(fX);

    fCstatu = "PROGRESS  ";

    if (fISW[4] >= 1) mnprin(5, fAmin);

    fEDM  = fBigedm;

    sig2  = fEDM;

    ncycl = 0;

//*-*-                                       . . . . .  start main loop

L50:

    if (sig2 < fEpsi && fEDM < fEpsi) goto L76;

    sig2 = fEDM;

    if (fNfcn - npfn > fNfcnmx) goto L78;

//*-*-        calculate new point * by reflection

    for (i = 1; i <= fNpar; ++i) {

	pb = 0;

	for (j = 1; j <= nparp1; ++j) { pb += wg*fP[i + j*fMaxpar - fMaxpar-1]; }

	fPbar[i-1]  = pb - wg*fP[i + jh*fMaxpar - fMaxpar-1];

	fPstar[i-1] = (alpha + 1)*fPbar[i-1] - alpha*fP[i + jh*fMaxpar - fMaxpar-1];

    }

    mninex(fPstar);

    (*fFCN)(nparx, fGin, ystar, fU, 4);    ++fNfcn;

    if (ystar >= fAmin) goto L70;

//*-*-        point * better than jl, calculate new point **

    for (i = 1; i <= fNpar; ++i) {

	fPstst[i-1] = gamma*fPstar[i-1] + (1 - gamma)*fPbar[i-1];

    }

    mninex(fPstst);

    (*fFCN)(nparx, fGin, ystst, fU, 4);    ++fNfcn;

//*-*-        try a parabola through ph, pstar, pstst.  min = prho

    y1 = (ystar - y[jh-1])*rho2;

    y2 = (ystst - y[jh-1])*rho1;

    rho = (rho2*y1 - rho1*y2)*.5 / (y1 - y2);

    if (rho < rhomin) goto L66;

    if (rho > rhomax) rho = rhomax;

    for (i = 1; i <= fNpar; ++i) {

	fPrho[i-1] = rho*fPbar[i-1] + (1 - rho)*fP[i + jh*fMaxpar - fMaxpar-1];

    }

    mninex(fPrho);

    (*fFCN)(nparx, fGin, yrho, fU, 4);    ++fNfcn;

    if (yrho < y[jl-1] && yrho < ystst) goto L65;

    if (ystst < y[jl-1]) goto L67;

    if (yrho > y[jl-1])  goto L66;

//*-*-        accept minimum point of parabola, PRHO

L65:

    mnrazz(yrho, fPrho, y, jh, jl);

    goto L68;

L66:

    if (ystst < y[jl-1]) goto L67;

    mnrazz(ystar, fPstar, y, jh, jl);

    goto L68;

L67:

    mnrazz(ystst, fPstst, y, jh, jl);

L68:

    ++ncycl;

    if (fISW[4] < 2) goto L50;

    if (fISW[4] >= 3 || ncycl % 10 == 0) {

	mnprin(5, fAmin);

    }

    goto L50;

//*-*-        point * is not as good as jl

L70:

    if (ystar >= y[jh-1]) goto L73;

    jhold = jh;

    mnrazz(ystar, fPstar, y, jh, jl);

    if (jhold != jh) goto L50;

//*-*-        calculate new point **

L73:

    for (i = 1; i <= fNpar; ++i) {

	fPstst[i-1] = beta*fP[i + jh*fMaxpar - fMaxpar-1] + (1 - beta)*fPbar[i-1];

    }

    mninex(fPstst);

    (*fFCN)(nparx, fGin, ystst, fU, 4);    ++fNfcn;

    if (ystst > y[jh-1]) goto L1;

//*-*-    point ** is better than jh

    if (ystst < fAmin) goto L67;

    mnrazz(ystst, fPstst, y, jh, jl);

    goto L50;

//*-*-                                       . . . . . .  end main loop

L76:

    if (fISW[4] >= 0) {

	fPrintf(" SIMPLEX MINIMIZATION HAS CONVERGED.");

    }

    fISW[3] = 1;

    goto L80;

L78:

    if (fISW[4] >= 0) {

	fPrintf(" SIMPLEX TERMINATES WITHOUT CONVERGENCE.");

    }

    fCstatu = "CALL LIMIT";

    fISW[3] = -1;

    fISW[0] = 1;

L80:

    for (i = 1; i <= fNpar; ++i) {

	pb = 0;

	for (j = 1; j <= nparp1; ++j) { pb += wg*fP[i + j*fMaxpar - fMaxpar-1]; }

	fPbar[i-1] = pb - wg*fP[i + jh*fMaxpar - fMaxpar-1];

    }

    mninex(fPbar);

    (*fFCN)(nparx, fGin, ypbar, fU, 4);    ++fNfcn;

    if (ypbar < fAmin) 	mnrazz(ypbar, fPbar, y, jh, jl);

    mninex(fX);

    if (fNfcnmx + npfn - fNfcn < fNpar*3) goto L90;

    if (fEDM > fEpsi*2) goto L1;

L90:

    if (fISW[4] >= 0) mnprin(5, fAmin);

} /* mnsimp_ */



//______________________________________________________________________________

void Midnight::mnstat(MDouble &fmin, MDouble &fedm, MDouble &errdef, MInt &npari, MInt &nparx, MInt &istat)

{



//*-*-*-*-*Returns concerning the current status of the minimization*-*-*-*-*

//*-*      =========================================================

//*-*       User-called

//*-*          Namely, it returns:

//*-*        FMIN: the best function value found so far

//*-*        FEDM: the estimated vertical distance remaining to minimum

//*-*        ERRDEF: the value of UP defining parameter uncertainties

//*-*        NPARI: the number of currently variable parameters

//*-*        NPARX: the highest (external) parameter number defined by user

//*-*        ISTAT: a status integer indicating how good is the covariance

//*-*           matrix:  0= not calculated at all

//*-*                    1= approximation only, not accurate

//*-*                    2= full matrix, but forced positive-definite

//*-*                    3= full accurate covariance matrix

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    fmin   = fAmin;

    fedm   = fEDM;

    errdef = fUp;

    npari  = fNpar;

    nparx  = fNu;

    istat  = fISW[1];

    if (fEDM == fBigedm) fedm = fUp;

    if (fAmin == fUndefi) {

	fmin  = 0;

	fedm  = fUp;

	istat = 0;

    }

} /* mnstat_ */



//______________________________________________________________________________

void Midnight::mntiny(MDouble epsp1, MDouble &epsbak)

{

//*-*-*-*-*-*-*-*To find the machine precision*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

//*-*            =============================

//*-*        Compares its argument with the value 1.0, and returns

//*-*        the value .TRUE. if they are equal.  To find EPSMAC

//*-*        safely by foiling the Fortran optimizer

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    epsbak = epsp1 - 1;

} /* mntiny_ */



//______________________________________________________________________________

MBool Midnight::mnunpt(MString &cfname)

{

//*-*-*-*-*-*Returns .TRUE. if CFNAME contains unprintable characters*-*-*-*

//*-*        ========================================================

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    static MInt i, l, ic;

    MBool ret_val;

    static MString cpt = " ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890./;:[]$%*_!@#&+()";



    ret_val = kFALSE;

    l       = strlen((const char*)cfname);

    for (i = 1; i <= l; ++i) {

	for (ic = 1; ic <= 80; ++ic) {

	    if (cfname[i-1] == cpt[ic-1]) goto L100;

	}

	return kTRUE;

L100:

	;

    }

    return ret_val;

} /* mnunpt_ */



//______________________________________________________________________________

void Midnight::mnvert(MDouble *a, MInt l, MInt, MInt n, MInt &ifail)

{

//*-*-*-*-*-*-*-*-*-*-*-*Inverts a symmetric matrix*-*-*-*-*-*-*-*-*-*-*-*-*

//*-*                    ==========================

//*-*        inverts a symmetric matrix.   matrix is first scaled to

//*-*        have all ones on the diagonal (equivalent to change of units)

//*-*        but no pivoting is done since matrix is positive-definite.

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    /* System generated locals */

    MInt a_offset;



    /* Local variables */

    static MDouble q[kMAXDIM], s[kMAXDIM], si, pp[kMAXDIM];

    static MInt i, j, k, kp1, km1;



    /* Parameter adjustments */

    a_offset = l + 1;

    a -= a_offset;



    /* Function Body */

    ifail = 0;

    if (n < 1) goto L100;

    if (n > fMaxint) goto L100;

//*-*-                  scale matrix by sqrt of diag elements

    for (i = 1; i <= n; ++i) {

	si = a[i + i*l];

	if (si <= 0) goto L100;

	s[i-1] = 1 / TMath_Sqrt(si);

    }

    for (i = 1; i <= n; ++i) {

	for (j = 1; j <= n; ++j) {

	    a[i + j*l] = a[i + j*l]*s[i-1]*s[j-1];

	}

    }

//*-*-                                       . . . start main loop . . . .

    for (i = 1; i <= n; ++i) {

	k = i;

//*-*-                  preparation for elimination step1

	if (a[k + k*l] != 0) q[k-1] = 1 / a[k + k*l];

	else goto L100;

	pp[k-1] = 1;

	a[k + k*l] = 0;

	kp1 = k + 1;

	km1 = k - 1;

	if (km1 < 0) goto L100;

	else if (km1 == 0) goto L50;

	else               goto L40;

L40:

	for (j = 1; j <= km1; ++j) {

	    pp[j-1]    = a[j + k*l];

	    q[j-1]     = a[j + k*l]*q[k-1];

	    a[j + k*l] = 0;

	}

L50:

	if (k - n < 0) goto L51;

	else if (k - n == 0) goto L60;

	else                goto L100;

L51:

	for (j = kp1; j <= n; ++j) {

	    pp[j-1]    = a[k + j*l];

	    q[j-1]     = -a[k + j*l]*q[k-1];

	    a[k + j*l] = 0;

	}

//*-*-                  elimination proper

L60:

	for (j = 1; j <= n; ++j) {

	    for (k = j; k <= n; ++k) { a[j + k*l] += pp[j-1]*q[k-1]; }

	}

    }

//*-*-                  elements of left diagonal and unscaling

    for (j = 1; j <= n; ++j) {

	for (k = 1; k <= j; ++k) {

	    a[k + j*l] = a[k + j*l]*s[k-1]*s[j-1];

	    a[j + k*l] = a[k + j*l];

	}

    }

    return;

//*-*-                  failure return

L100:

    ifail = 1;

} /* mnvert_ */



//______________________________________________________________________________

void Midnight::mnwarn(const char *copt1, const char *corg1, const char *cmes1)

{

//*-*-*-*-*-*-*-*-*-*-*-*Prints Warning messages*-*-*-*-*-*-*-*-*-*-*-*-*-*

//*-*                    =======================

//*-*     If COPT='W', CMES is a WARning message from CORG.

//*-*     If COPT='D', CMES is a DEBug message from CORG.

//*-*         If SET WARnings is in effect (the default), this routine

//*-*             prints the warning message CMES coming from CORG.

//*-*         If SET NOWarnings is in effect, the warning message is

//*-*             stored in a circular buffer of length MAXMES.

//*-*         If called with CORG=CMES='SHO', it prints the messages in

//*-*             the circular buffer, FIFO, and empties the buffer.

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    MString copt = copt1;

    MString corg = corg1;

    MString cmes = cmes1;



    const MInt MAXMES = 10;

    static MInt ityp, i, ic, nm;

    static MString englsh, ctyp;



    if (corg(0,3) != "SHO" || cmes(0,3) != "SHO") {



//*-*-            Either print warning or put in buffer

       if (copt == "W") {

          ityp = 1;

	  if (fLwarn) {

             fPrintf(" MINUIT WARNING IN %s",(const char*)corg);

             fPrintf(" ============== ",(const char*)cmes);

	     return;

	  }

       } else {

	  ityp = 2;

	  if (fLrepor) {

             fPrintf(" MINUIT DEBUG FOR %s",(const char*)corg);

             fPrintf(" =============== %s ",(const char*)cmes);

	     return;

	  }

       }

//*-*-                if appropriate flag is off, fill circular buffer

       if (fNwrmes[ityp-1] == 0) fIcirc[ityp-1] = 0;

       ++fNwrmes[ityp-1];

       ++fIcirc[ityp-1];

       if (fIcirc[ityp-1] > 10) 	fIcirc[ityp-1] = 1;

       ic = fIcirc[ityp-1];

       fOrigin[ic + ityp*10 - 11] = corg;

       fWarmes[ic + ityp*10 - 11] = cmes;

       fNfcwar[ic + ityp*10 - 11] = fNfcn;

       return;

   }



//*-*-            'SHO WARnings', ask if any suppressed mess in buffer

    if (copt == "W") {

	ityp = 1;

	ctyp = "WARNING";

    } else {

	ityp = 2;

	ctyp = "*DEBUG*";

    }

    if (fNwrmes[ityp-1] > 0) {

	englsh = " WAS SUPPRESSED.  ";

	if (fNwrmes[ityp-1] > 1) englsh = "S WERE SUPPRESSED.";

        fPrintf(" %5d MINUIT %s MESSAGE%s",fNwrmes[ityp-1]

             ,(const char*)ctyp,(const char*)englsh);

	nm = fNwrmes[ityp-1];

	ic = 0;

	if (nm > MAXMES) {

	    fPrintf(" ONLY THE MOST RECENT 10 WILL BE LISTED BELOW.");

	    nm = MAXMES;

	    ic = fIcirc[ityp-1];

	}

        fPrintf("  CALLS  ORIGIN         MESSAGE");

	for (i = 1; i <= nm; ++i) {

	    ++ic;

	    if (ic > MAXMES) ic = 1;

            fPrintf(" %6d %s %s", fNfcwar[ic + ityp*10 - 11],

                               (const char*)fOrigin + (ic + ityp*10 - 11)*10,

                               (const char*)fWarmes + (ic + ityp*10 - 11)*60);

	}

	fNwrmes[ityp-1] = 0;

	fPrintf(" ");

    }

} /* mnwarn_ */



//______________________________________________________________________________

void Midnight::mnwerr()

{

//*-*-*-*-*-*-*-*Calculates the WERR, external parameter errors*-*-*-*-*-*-*

//*-*            ==============================================

//*-*      and the global correlation coefficients, to be called

//*-*      whenever a new covariance matrix is available.

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*



    static MDouble denom, ba, al, dx, du1, du2;

    static MInt ndex, ierr, i, j, k, l, ndiag, k1, iin;



//*-*-                        calculate external error if v exists

    if (fISW[1] >= 1) {

	for (l = 1; l <= fNpar; ++l) {

	    ndex = l*(l + 1) / 2;

	    dx = TMath_Sqrt(TMath_Abs(fVhmat[ndex-1]*fUp));

	    i = fNexofi[l-1];

	    if (fNvarl[i-1] > 1) {

		al = fAlim[i-1];

		ba = fBlim[i-1] - al;

		du1 = al + 0.5*(TMath_Sin(fX[l-1] + dx) + 1)*ba - fU[i-1];

		du2 = al + 0.5*(TMath_Sin(fX[l-1] - dx) + 1)*ba - fU[i-1];

		if (dx > 1) du1 = ba;

		dx = 0.5*(TMath_Abs(du1) + TMath_Abs(du2));

	    }

	    fWerr[l-1] = dx;

	}

    }

//*-*-                         global correlation coefficients

    if (fISW[1] >= 1) {

	for (i = 1; i <= fNpar; ++i) {

	    fGlobcc[i-1] = 0;

	    k1 = i*(i-1) / 2;

	    for (j = 1; j <= i; ++j) {

		k = k1 + j;

		fP[i + j*fMaxpar - fMaxpar-1] = fVhmat[k-1];

		fP[j + i*fMaxpar - fMaxpar-1] = fP[i + j*fMaxpar - fMaxpar-1];

	    }

	}

	mnvert(fP, fMaxint, fMaxint, fNpar, ierr);

	if (ierr == 0) {

	    for (iin = 1; iin <= fNpar; ++iin) {

		ndiag = iin*(iin + 1) / 2;

		denom = fP[iin + iin*fMaxpar - fMaxpar-1]*fVhmat[ndiag-1];

		if (denom <= 1 && denom >= 0) fGlobcc[iin-1] = 0;

		else                          fGlobcc[iin-1] = TMath_Sqrt(1 - 1 / denom);

	    }

	}

    }

} /* mnwerr_ */

//////////////////////////////////////////////////////////////////////////////

static char* Form(

 char* aFormat

,...

)

//////////////////////////////////////////////////////////////////////////////

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//

{

  va_list args;

  va_start(args,aFormat);

  vsprintf(sBuffer,aFormat,args);

  va_end(args);

  return sBuffer;

}



//////////////////////////////////////////////////////////////////////////////

static void Printf(

 const char* aFormat

,...

)

//////////////////////////////////////////////////////////////////////////////

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//

{

  va_list args;

  va_start(args,aFormat);

  vprintf(aFormat,args);

  va_end(args);

  printf("\n");

}

