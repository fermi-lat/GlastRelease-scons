
//
//   Example of a program to fit non-equidistant data points
//   =======================================================
//
//   The fitting function fcn is a simple chisquare function
//   The data consists of 4 data points (arrays x,y) + the errors in errorsy
//   More details on the various functions or parameters for these functions
//   can be obtained in an interactive ROOT session with:
//    Root > TMinuit *minuit = new TMinuit(10);
//    Root > minuit->mnhelp("*")  to see the list of possible keywords
//    Root > minuit->mnhelp("SET") explains most parameters
//

Float_t x[4],y[4],errory[4];
Double_t param[2];

//________
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
   const Int_t nbins = 4;
   Int_t i;

//calculate chisquare
   Double_t chisq = 0;
   Double_t delta;
   for (i=0;i<nbins; i++) {
     delta  = (y[i]-func(x[i],par))/errory[i];
     chisq += delta*delta;
   }
   f = chisq;
}

//_______ Simple linear function
Double_t func(Float_t x,Double_t *par)
{
 Double_t value=(x*par[0] + par[1]);
 return value;
}

//________
void DoFit()
{
// The y values
        y[0]=1.3;
        y[1]=0.47;
        y[2]=0.01;
        y[3]=-0.28;
// The errors on y values
        errory[0]=0.1;
        errory[1]=0.2;
        errory[2]=0.6;
        errory[3]=0.1;
// the x values
        x[0]=1.0;
        x[1]=2.0;
        x[2]=3.0;
        x[3]=4.0;


   TMinuit *gMinuit = new TMinuit(5);  //initialize TMinuit with a maximum of 5 params
   gMinuit->SetFCN(fcn);

   Double_t arglist[10];
   Int_t ierflg = 0;

   arglist[0] = 1;
   gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

// Set starting values and step sizes for parameters
   static Double_t vstart[2] = {-0.5, 1.0 };
   static Double_t step[2] = {0.1 , 0.1};
   gMinuit->mnparm(0, "a1", vstart[0], step[0], 0,0,ierflg);
   gMinuit->mnparm(1, "a2", vstart[1], step[1], 0,0,ierflg);


// Now ready for minimization step
   arglist[0] = 500;
   arglist[1] = 1.;
   gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

// Print results
   Double_t amin,edm,errdef;
   Int_t nvpar,nparx,icstat;
   gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
   gMinuit->mnprin(3,amin);

// Save the fitted parameters

   Double_t temp;
   gMinuit->GetParameter(0,param[0],temp);
   gMinuit->GetParameter(1,param[1],temp);

}

void DisplayFit(){
   c1 = new TCanvas("c1","Minuit Linear Fit ",200,10,700,500);
   c1->SetGridx();
   c1->SetGridy();

   // Draw the points
   TGraphErrors * g = new TGraphErrors(4,x,y,0,errory);
   g->SetMarkerColor(4);
   g->SetMarkerStyle(21);
   g->Draw("ALP");

   // Create the function and draw it
   fun1 = new TF1("fun1","x*[0] + [1]",0.0,5.0);
   fun1->SetParameters(param);
   fun1->Draw("SAME");

   c1->Update();

}
