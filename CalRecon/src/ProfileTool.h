
#ifndef __ProfileTool_H
#define __ProfileTool_H 1

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "EnergyCorr.h"
//Gamma function and Minuit
#include "TMath.h"
#include "TMinuit.h"

/**   
* @class ProfileTool
*
* Algorithm for calculating energy leakage by fitting the longitudinal
* shower profile.
*
*
* $Header$
*/


class ProfileTool :  public EnergyCorr {

public:
    
    //! destructor
    ProfileTool( const std::string& type, const std::string& name, const IInterface* parent);
     ~ProfileTool() {}; 
    
     StatusCode initialize();

        

//! Longitudinal profile fitting method
/*! It performs a longitudinal profile fitting using the incomplete gamma
* function ( gamma.cxx). The mean energy density per length unit is taken as:
* \f[
f(z) = \frac{1}{\lambda}
\frac{exp(-z/\lambda)(z/\lambda)^{(1-\alpha)}}{\Gamma(\alpha)}
\f]
* To take into account the angle \f$ \theta \f$ of the particle trajectory
* with vertical direction, z should be divided by \f$ cos(\theta) \f$ . 
  
* Thus integrated on the i-th crystal pathlength it gives 
* \f[
E_i = E_{tot}(\Gamma_{inc}
(z_i/\lambda/cos(\theta),\alpha)-
\Gamma_{inc}(z_{i+1}/\lambda/cos(\theta),\alpha))
\f]
* The 2 shower parameters here \f$ \alpha \f$ and \f$ \lambda \f$ describe the maximum
* position and the exponential decrease of the profile. Those parameters
* have been estimated using
* a MC and their dependance over E has been fitted by a power law.
* They are log-normally 
* distributed with a very broad distribution
* ( hence some of the shower fluctuations ).
* Therefore they should not be included in the fitting process. 
*
* Here we use 4 parameters: total energy,
*                           starting point,
*                           \f$ \alpha \f$ and
*                           \f$ \lambda \f$.
*
* The Minuit package, used for fitting, permits to fix any parameter
* or to include it into fitting process. In the current implementation
* of Profile() function \f$ \alpha \f$ and \f$ \lambda \f$ are fixed
* at the values calculated by predfined formulas from uncorrected total energy
* sum in the calorimeter, while 2 other parameters are optimized by Minuit
* to get the best fit.
*
* The fit output energy is debiased using a polynomial:
*     P( log( fit energy ), cos( angle of incdence ) ) giving the most probable
* bias for a certain monte carlo energy and incidfent angle. Its parameters are 
* contained in \$(CALRECONROOT)/xml/CalProfile.xml
* The correction is operated thrice:
*   \f[ E_{debiased}=E_{fit} 
    \f]
    \f[ E_{debiased}=E_{fit}+E_{debiased}*\mathrm{P}(E_{debiased},\cos(\theta))
        \;\rbrace \times 3 
    \f] 
*
* The input is:
* \param eTotal total energy measured in the calorimeter in MeV
* \param cl     the CalCluster in which the results are saved
*
* The output ( the 4 fitting parameters and the chi square) of this method is
* stored in the CalCluster.  
* 
* \warning Gives sensible results only at large enough energies ie ~10GeV
*
* \b Revision:
* - 10/17/00    RT    comments added
* - 05/00       RT    first implementation
*/     
     StatusCode doEnergyCorr(double eTotal, Event::CalCluster* cluster);

     StatusCode finalize();
    
     StatusCode execute();

protected:

     // function passed to Minuit to minimize
    static void fcn(int & , double *, double &f, double *par, int );
    // gamma function for shower shape 
    static double gam_prof(double *par, int i);
    // profile bias function
    double bias( double logEnergy );

private:

    TMinuit* m_minuit;

    static double m_xtalHeight; //!< xtal height in cm
    static double m_xtalWidth;  //!< xtal width  in cm
    static int m_nbins;  //!< Number of bins used for the fit
    static std::vector<double> m_g_elayer;  //!< Energy per layer in GeV

    /// input XML file containing parameters for bias correction
    std::string	m_xmlFile;
    std::vector<double> m_BiasVector;
    int m_BiasParNb;
    int m_BiasCTNb;
    double m_BiasCTLim;
};

#endif



