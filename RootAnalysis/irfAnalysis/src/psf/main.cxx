/** @file main.h
  $Header$
  */
  
#include "PSF.h"
#include "MakeDists.h"
#include "SumOfGaussians.h"
#include "LinearModel.h"
#include "PsfModel.h"

namespace {
   void setBounds(Fitter * fitter) {
// Set bounds on the sum-of-Gaussians fit parameters.
      double lower[] = {0., -3.5, 0., 0., -3.5, 0.};
      double upper[] = {1e3, -0.5, 3., 1e3, -0.5, 3.};
      std::vector<double> lower_bounds(lower, lower+6);
      std::vector<double> upper_bounds(upper, upper+6);
      fitter->setBounds(lower_bounds, upper_bounds);
   }
}

int main(){

    PSF p("psf.root");
    if(!p.fileExists())  p.project( 0, 5.0, 50);
    std::string ps = "psf.ps";
    p.draw(ps+"(", 0.2);
    p.drawAeff(ps);
    p.drawError(ps);
    p.drawAsymmetry(ps+")");

// Fit the angle ErrDists with a sum of two Gaussian functions.
    Fitter * twoGauss = new SumOfGaussians("Tkr1ThetaErrFits.root");
    ::setBounds(twoGauss);

// Theta error.
    MakeDists thetaErrDist("Tkr1ThetaErr.root");
    if(! thetaErrDist.fileExists() )
        thetaErrDist.project("log10(Tkr1ThetaErr)", -3.5, -0.5, 50, twoGauss);
    thetaErrDist.draw( "Tkr1ThetaErr.ps", 0.3 );
    delete twoGauss;
    thetaErrDist.addCutInfo("Tkr1ThetaErrFits.root", "fitParams");

// Phi error.
    twoGauss = new SumOfGaussians("Tkr1PhiErrFits.root");
    ::setBounds(twoGauss);
    MakeDists phiErrDist("Tkr1PhiErr.root");
    if(! phiErrDist.fileExists() )
        phiErrDist.project("log10(Tkr1PhiErr)", -3.5, -0.5, 50, twoGauss);
    phiErrDist.draw("Tkr1PhiErr.ps", 0.3);
    delete twoGauss;
    phiErrDist.addCutInfo("Tkr1PhiErrFits.root", "fitParams");

// Theta and phi errors added in quadrature.
    twoGauss = new SumOfGaussians("Tkr1ErrFits.root");
    ::setBounds(twoGauss);
    MakeDists errDist("Tkr1Err.root");
    if(! errDist.fileExists() )
       errDist.project("log10(sqrt(pow(Tkr1ThetaErr,2.)+pow(Tkr1PhiErr,2.)))",
                       -3.5, 0., 50, twoGauss);
    errDist.draw("Tkr1Err.ps", 0.3);
    delete twoGauss;
    errDist.addCutInfo("Tkr1ErrFits.root", "fitParams");

// Fit the profile of log10(Tkr1PhiErr) vs log10(Tkr1ThetaErr).
    Fitter * linearModel = new LinearModel("Tkr1ErrProfileFits.root");
    bool makeProfile(true);
    MakeDists tkr1ErrProfile("Tkr1ErrProfile.root", makeProfile);
    if(! tkr1ErrProfile.fileExists() )
       tkr1ErrProfile.project("log10(Tkr1ThetaErr):log10(Tkr1PhiErr)",
                              -3.5, -0.5, 100, linearModel);
    tkr1ErrProfile.draw("Tkr1ErrProfiles.ps", -0.5);
    delete linearModel;
    tkr1ErrProfile.addCutInfo("Tkr1ErrProfileFits.root", "fitParams");

// Fit the scaled PSF distributions using the function proposed by
// Luis Reyes and Steve Ritz.
    Fitter * psfModel = new PsfModel("scaledPsfFits.root");
    MakeDists scaledPsf("scaledPsf.root");
    if ( !scaledPsf.fileExists() )
       scaledPsf.project("BestDirErr/PSFscaleFactor", 0, 5, 100, psfModel);
    scaledPsf.draw("scaledPsf.ps", 0.1);
    delete psfModel;
    scaledPsf.addCutInfo("scaledPsfFits.root", "fitParams");

    std::cout << "done" << std::endl;
    return 0;
}

