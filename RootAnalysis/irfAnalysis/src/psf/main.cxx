/** @file main.h
  $Header$
  */
  
#include "PSF.h"
#include "MakeDists.h"
#include "SumOfGaussians.h"

int main(){

    PSF p;
    if(!p.fileExists())  p.project();
    std::string ps = "psf.ps";
    p.draw(ps+"(");
    p.drawAeff(ps);
    p.drawError(ps);
    p.drawAsymmetry(ps+")");
// Fit the angle ErrDists with a sum of two Gaussian functions.
    Fitter * twoGauss = new SumOfGaussians();

// Set bounds on fit parameters.
    double lower[] = {0., -3.5, 0., 0., -3.5, 0.};
    double upper[] = {1e3, -0.5, 3., 1e3, -0.5, 3.};
    std::vector<double> lower_bounds(lower, lower+6);
    std::vector<double> upper_bounds(upper, upper+6);
    twoGauss->setBounds(lower_bounds, upper_bounds);

// Create the distributions.
    MakeDists thetaErrDist("Tkr1ThetaErr.root");
    thetaErrDist.project("log10(Tkr1ThetaErr)", -3.5, -0.5, 50, twoGauss);
    thetaErrDist.draw( "Tkr1ThetaErr.ps", 0.3 );

    MakeDists phiErrDist("Tkr1PhiErr.root");
    phiErrDist.project("log10(Tkr1PhiErr)", -3.5, -0.5, 50, twoGauss);
    phiErrDist.draw("Tkr1PhiErr.ps", 0.3);

    std::cout << "done" << std::endl;
    return 0;
}
