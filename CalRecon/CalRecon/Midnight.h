//  The below code is a modified version of the ROOT/TMinuit class
// in order to have a stand alone C++ version of the Minuit package.
//    G.Barrand

#ifndef Midnight_h
#define Midnight_h

#include "CalRecon/MString.h"

typedef unsigned char MBool;
typedef int MInt;
typedef double MDouble;
typedef void(*MFunction)(MInt&,MDouble*,MDouble&,MDouble*,MInt);
typedef void(*MPrintf)(const char*,...);

class Midnight {
private:
        MInt        fEmpty;            //Initialization flag (1 = Midnight initialized)
        MInt        fMaxpar;           //Maximum number of parameters
        MString      *fCpnam;           //Array of parameters names
        MDouble     *fU;               //External (visible to user in FCN) value of parameters
        MDouble     *fAlim;            //Lower limits for parameters. If zero no limits
        MDouble     *fBlim;            //Upper limits for parameters
        MDouble     *fErp;             //Positive Minos errors if calculated
        MDouble     *fErn;             //Negative Minos errors if calculated
        MDouble     *fWerr;            //External parameters error (standard deviation, defined by UP)
        MDouble     *fGlobcc;          //Global Correlation Coefficients
           MInt     *fNvarl;           //parameters flag (-1=undefined, 0=constant..)
           MInt     *fNiofex;          //Internal parameters number, or zero if not currently variable
           MInt     *fNexofi;          //External parameters number for currently variable parameters
        MDouble     *fX;               //Internal parameters values
        MDouble     *fXt;              //Internal parameters values X saved as Xt
        MDouble     *fDirin;           //(Internal) step sizes for current step
        MDouble     *fXs;              //Internal parameters values saved for fixed params
        MDouble     *fXts;             //Internal parameters values X saved as Xt for fixed params
        MDouble     *fDirins;          //(Internal) step sizes for current step for fixed params
        MDouble     *fGrd;             //First derivatives
        MDouble     *fG2;              //
        MDouble     *fGstep;           //Step sizes
        MDouble     *fGin;             //
        MDouble     *fDgrd;            //Uncertainties
        MDouble     *fGrds;            //
        MDouble     *fG2s;             //
        MDouble     *fGsteps;          //
        MInt        *fIpfix;           //List of fixed parameters
        MInt        fNpfix;            //Number of fixed parameters
        MDouble     *fVhmat;           //(Internal) error matrix stored as Half MATrix, since it is symmetric
        MDouble     *fVthmat;          //VHMAT is sometimes saved in VTHMAT, especially in MNMNOT
        MDouble     *fP;               //
        MDouble     *fPstar;           //
        MDouble     *fPstst;           //
        MDouble     *fPbar;            //
        MDouble     *fPrho;            //Minimum point of parabola
        MInt        fMaxint;           //Maximum number of internal parameters
        MInt        fNpar;             //Number of parameters
        MInt        fMaxext;           //Maximum number of external parameters
        MInt        fNu;               //
        MInt        fIsysrd;           //standardInput unit
        MInt        fIsyswr;           //standard output unit
        MInt        fIsyssa;           //
        MInt        fNpagwd;           //Page width
        MInt        fNpagln;           //Number of lines per page
        MInt        fNewpag;           //
        MInt        fIstkrd[10];       //
        MInt        fNstkrd;           //
        MInt        fIstkwr[10];       //
        MInt        fNstkwr;           //
        MString      fCfrom;            //
        MString      fCstatu;           //
        MString      fCtitl;            //
        MString      fCword;            //
        MString      fCundef;           //
        MString      fCvrsn;            //
        MString      fCovmes[4];        //
        MInt        fISW[7];           //Array of switches
        MInt        fIdbg[11];         //Array of internal debug switches
        MInt        fNblock;           //Number of Minuit data blocks
        MInt        fIcomnd;           //Number of commands
        MDouble     fAmin;             //Minimum value found for FCN
        MDouble     fUp;               //FCN+-UP defines errors (for chisquare fits UP=1)
        MDouble     fEDM;              //Estimated vertical distance to the minimum
        MDouble     fFval3;            //
        MDouble     fEpsi;             //
        MDouble     fApsi;             //
        MDouble     fDcovar;           //Relative change in covariance matrix
        MInt        fNfcn;             //Number of calls to FCN
        MInt        fNfcnmx;           //Maximum number of calls to FCN
        MInt        fNfcnlc;           //
        MInt        fNfcnfr;           //
        MInt        fItaur;            //
        MInt        fIstrat;           //
        MInt        fNwrmes[2];        //
        MDouble     *fWord7;           //
        MBool       fLwarn;            //true if warning messges are to be put out (default=true)
        MBool       fLrepor;           //true if exceptional conditions are put out (default=false)
        MBool       fLimset;           //true if a parameter is up against limits (for MINOS)
        MBool       fLnolim;           //true if there are no limits on any parameters (not yet used)
        MBool       fLnewmn;           //true if the previous process has unexpectedly improved FCN
        MBool       fLphead;           //true if a heading should be put out for the next parameter definition
        MDouble     fEpsmac;           //machine precision for floating points:
        MDouble     fEpsma2;           //sqrt(fEpsmac)
        MDouble     fVlimlo;           //
        MDouble     fVlimhi;           //
        MDouble     fUndefi;           //Undefined number = -54321
        MDouble     fBigedm;           //Big EDM = 123456
        MDouble     fUpdflt;           //
        MDouble     *fXpt;             //X array of points for contours
        MDouble     *fYpt;             //Y array of points for contours
        MString      *fChpt;            //Character to be plotted at the X,Y contour positions
        MDouble     fXmidcr;           //
        MDouble     fYmidcr;           //
        MDouble     fXdircr;           //
        MDouble     fYdircr;           //
        MInt        fKe1cr;            //
        MInt        fKe2cr;            //
        MString      *fOrigin;          //
        MString      *fWarmes;          //
        MInt        fNfcwar[20];       //
        MInt        fIcirc[2];         //
        MFunction fFCN;
        MPrintf fPrintf;


// methods performed on Midnight class
private:
	   Midnight(const Midnight &m);
 void      BuildArrays(MInt maxpar=15);
 void      DeleteArrays();
public:
           Midnight();
           Midnight(MInt maxpar);
 virtual   ~Midnight();
 void      SetFCN(MFunction);
 void      SetPrintf(MPrintf);
 //

 // added R Terrier
 int GetParameter(int parNo, double &currentValue, double &currentError ); 


 void      mnamin();
 void      mnbins(MDouble a1, MDouble a2, MInt naa, MDouble &bl, MDouble &bh, MInt &nb, MDouble &bwid);
 void      mncalf(MDouble *pvec, MDouble &ycalf);
 void      mncler();
 void      mncntr(MInt ke1, MInt ke2, MInt &ierrf);
 void      mncomd(MString crdbin, MInt &icondn);
 void      mncont(MInt ke1, MInt ke2, MInt nptu, MDouble *xptu, MDouble *yptu, MInt &ierrf);
 void      mncrck(MString &crdbuf, MInt maxcwd, MString &comand, MInt &lnc
                ,  MInt mxp, MDouble *plist, MInt &llist, MInt &ierr, MInt isyswr);
 void      mncros(MDouble &aopt, MInt &iercr);
 void      mncuve();
 void      mnderi();
 void      mndxdi(MDouble pint, MInt ipar, MDouble &dxdi);
 void      mneig(MDouble *a, MInt ndima, MInt n, MInt mits, MDouble *work, MDouble precis, MInt &ifault);
 void      mnemat(MDouble *emat, MInt ndim);
 void      mnerrs(MInt number, MDouble &eplus, MDouble &eminus, MDouble &eparab, MDouble &gcc);
 void      mneval(MDouble anext, MDouble &fnext, MInt &ierev);
 void      mnexcm(MString comand, MDouble *plist, MInt llist, MInt &ierflg) ;
 void      mnexin(MDouble *pint);
 void      mnfixp(MInt iint, MInt &ierr);
 void      mnfree(MInt k);
 void      mngrad();
 void      mnhelp(MString comd);
 void      mnhess();
 void      mnhes1();
 void      mnimpr();
 void      mninex(MDouble *pint);
 void      mninit(MInt i1, MInt i2, MInt i3);
 void      mnlims();
 void      mnline(MDouble *start, MDouble fstart, MDouble *step, MDouble slope, MDouble toler);
 void      mnmatu(MInt kode);
 void      mnmigr();
 void      mnmnos();
 void      mnmnot(MInt ilax, MInt ilax2, MDouble &val2pl, MDouble &val2mi);
 void      mnparm(MInt k, MString cnamj, MDouble uk, MDouble wk, MDouble a, MDouble b, MInt &ierflg);
 void      mnpars(MString &crdbuf, MInt &icondn);
 void      mnpfit(MDouble *parx2p, MDouble *pary2p, MInt npar2p, MDouble *coef2p, MDouble &sdev2p);
 void      mnpint(MDouble &pexti, MInt i, MDouble &pinti);
 void      mnplot(MDouble *xpt, MDouble *ypt, MString *chpt, MInt nxypt, MInt npagwd, MInt npagln);
 void      mnpout(MInt iuext, MString& chnam, MDouble &val, MDouble &err, MDouble &xlolim, MDouble &xuplim, MInt &iuint);
 void      mnprin(MInt inkode, MDouble fval);
 void      mnpsdf();
 void      mnrazz(MDouble ynew, MDouble *pnew, MDouble *y, MInt &jh, MInt &jl);
 void      mnrn15(MDouble &val, MInt &inseed);
 void      mnrset(MInt iopt);
 void      mnsave();
 void      mnscan();
 void      mnseek();
 void      mnset();
 void      mnsimp();
 void      mnstat(MDouble &fmin, MDouble &fedm, MDouble &errdef, MInt &npari, MInt &nparx, MInt &istat);
 void      mntiny(MDouble epsp1, MDouble &epsbak);
 MBool    mnunpt(MString &cfname);
 void      mnvert(MDouble *a, MInt l, MInt m, MInt n, MInt &ifail);
 void      mnwarn(const char *copt, const char *corg, const char *cmes);
 void      mnwerr();
};

#endif








