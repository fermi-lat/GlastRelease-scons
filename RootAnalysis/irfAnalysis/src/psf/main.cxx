/** @file main.h
  $Header$
  */
  
#include "PSF.h"
#include "MakeDists.h"
#include "SumOfGaussians.h"

namespace {
   void setBounds(Fitter * fitter) {
// Set bounds on fit parameters.
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

// Create the distributions.
    MakeDists thetaErrDist("Tkr1ThetaErr.root");
    if(! thetaErrDist.fileExists() )
        thetaErrDist.project("log10(Tkr1ThetaErr)", -3.5, -0.5, 50, twoGauss);
    thetaErrDist.draw( "Tkr1ThetaErr.ps", 0.3 );
    delete twoGauss;

    twoGauss = new SumOfGaussians("Tkr1PhiErrFits.root");
    ::setBounds(twoGauss);

    MakeDists phiErrDist("Tkr1PhiErr.root");
    if(! phiErrDist.fileExists() )
        phiErrDist.project("log10(Tkr1PhiErr)", -3.5, -0.5, 50, twoGauss);
    phiErrDist.draw("Tkr1PhiErr.ps", 0.3);
    delete twoGauss;

    std::cout << "done" << std::endl;
    return 0;
}
