
#ifndef __CALCLUSTERSALG_H
#define __CALCLUSTERSALG_H 1

#include "GaudiKernel/Algorithm.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "geometry/Vector.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/CalCluster.h"
class Midnight;

/**   
* @class CalClustersAlg
*
* Algorithm for reconstruction of energy and direction of incident particle
*
*
* Performs high level energy corrections
*
* The reconstruction here uses CalXtalRecCol to produce a CalClusterCol.
*  It evaluates the barycenter for each layer using the coordinates stored in
*  the CalXtalRecCol, and tries to correct for energy leakage using two
*  different methods:
*  -  Profile performs profile fitting
*  -  Leak performs a correction using the last layer correlation
* 
* See the description of Profile() and Leak() for details
* For most of the test beam energies Leak() should give better results.
* there can be a strong bias though because of miscalibration
*
* For a comparison one can see on the following plots the results of this
* method on R138 data and a MC run of 20 GeV positrons
* 
* \image html figurearticle2.gif
* \image latex figurearticle2.eps width=10cm
* \image html figuresim2.gif
* \image latex figuresim2.eps width=10cm
* \author Regis Terrier
* \author Alexandre Chekhtman
* \author Jose Hernando
* 
* \bug These energy correction methods are implemented on axis.
* \warning May not give sensible results when there is a large uncertainty in 
* the gains especially for the Leak method if there is a problem in the last 
* layer.
* \warning High energy corrections are intended for high energy 
* \todo Add low energy corrections 
*
*
* $Header$
*/


class CalClustersAlg : public Algorithm
{
public:
    
    //! constructor
    CalClustersAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    //! destructor
    virtual ~CalClustersAlg() {}; 
    
    StatusCode initialize();

        
/*!Performs the reconstruction. 
 * - Energy per layer is computed and stored in CalCluster in MeV
 * - Barycenter per layer is also computed and stored in CsICluster
 * - Performs high energy corrections: see Profile() and Leak() for details
 */        
    StatusCode execute();

/*! Finalization of algorithm 
 *  - minuit deleted
 */    
    StatusCode finalize();
    

    
//! Leakage correction using last layer correlation
/*!The second method uses the correlation between the escaping energy and the energy
*  deposited in the last layer of the calorimeter.  Indeed, the last layer carries
* the most important information concerning the leaking energy:  the total number
* of particles escaping through the back should be nearly proportional to the 
* energy deposited in the last layer.  The measured signal in that layer can therefore
* be modified to account for the leaking energy. 
*
* We used the Monte Carlo simulation of the GLAST beam test configuration to deter
* mine this correlation at several energies, from 2 GeV up to 40 GeV.  For one par
* ticular incident energy, the bidimensionnal distribution of the energy escaping 
* and the energy deposited in the last layer can be fitted by a simple linear function:
* \f[ E_{leak} = \alpha E_{last} + \beta \f]
* The \f$\alpha\f$ and \f$\beta\f$ parameters are proportional to the logarithm of the 
* incident energy and to its square, respectively.  Because the only information we 
* have, initially, on the incident energy is the total measured energy \f$E_m\f$, we have
* to use it as the estimator of \f$E_0\f$.  The reconstructed energy is then:
* \f[ E_{rec} = E_m + \alpha(E_m) E_{last} + \beta(E_m) \f]
* To improve the result, one can iterate using the new estimator to determine the 
* correct values of \f$\alpha\f$ and \f$\beta\f$.  
*
*\par The methods takes 2 arguments:
*\param eTotal Total energy measured in the calorimeter in MeV
*\param elast  Total energy measured in the last layer in MeV
*
*\return Corrected energy in MeV
*
* \warning Has been developped on axis only
* \warning Give biased results when the last layergains are misaligned
*
* \author Regis Terrier
*
* \b Revision: 
* - 10/17/00    RT    comments added
* - 09/15/00    RT    bias correction added, tuned on MC data
* - 08/20/00    RT    first implementation
*/
           
    double Leak(double sum,double elast);


//! Longitudinal profile fitting method
/*! It performs a longitudinal profile fitting using the incomplete gamma function
* ( gamma.cxx). The mean energy density per length unit is taken as:
* \f[
f(z) = \frac{1}{\lambda} \frac{exp(-z/\lambda)(z/\lambda)^{(1-\alpha)}}{\Gamma(\alpha)}
\f]
* Thus integrated on the i th Xtal pathlength it gives 
* \f[
E_i = E_{tot}(\Gamma_{inc}(z_i/\lambda,\alpha)-\Gamma_{inc}(z_{i+1}/\lambda,\alpha))
\f]
* the 2 shower parameters here alpha and lambda describes the maximum position and 
* and the exponential decrease of the profile. Those parameters have been estimated using
* a MC and their dependance over E has been fitted by a power law. They are log-normally 
* distributed with a very broad distribution ( hence some of the shower fluctuations ).
* Therefore they should not be included in the fitting process. 
*
* Here we use 4 parameters: total energy
*                           starting point
*                           alpha
*                           lambda
* They can be fixed or released in Profile()
*
*
* The input is:
* \param eTotal total energy measured in the calorimeter in MeV
* \param cl     the CsICluster in which the results are saved
*
* The output ( the 4 fitting parameters and the chi square) of this method is
* stored in the CsICluster.  
* 
* \warning The algorithm works only on axis thats why slope=1 in execute()
* \warning Gives sensible results only at large enough energies ie ~10GeV
*
* \b Revision:
* - 10/17/00    RT    comments added
* - 05/00       RT    first implementation
*/
    void Profile(double sum, Event::CalCluster* cl);


//! Reconstruct the direction in the calorimeter
/*! Basic algorithm for now, since we need to have knowledge
 *  on longitudinal errors
 *  Simply reconstruct direction on both sides XZ and YZ
 */
    Vector Fit_Direction(std::vector<Vector> pos,std::vector<Vector> sigma2,int nlayers);
    
protected:
    
    StatusCode retrieve();
    
private:
    
    //! reconstructed data for crystals, the input of the reconstruction
    Event::CalXtalRecCol* m_calXtalRecCol;
    //! the clusters list, the output of the reconstruction
    Event::CalClusterCol* m_calClusterCol;
    //! the minimizer for Profile()
    Midnight* minuit;

    //! this parameter permits to distinguish multiple calls
    //! to calorimeter reconstruction for the same event
    int m_callNumber;

    //! crystal width
    double m_CsIWidth;

    //! crystal height
    double m_CsIHeight;

    //! number of layers
    int m_CalnLayers;

    //! pointer to GlasDetSvc
    IGlastDetSvc* detSvc;

};

#endif



