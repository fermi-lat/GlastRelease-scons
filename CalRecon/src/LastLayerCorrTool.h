
#ifndef __LastLayerCorrTool_H
#define __LastLayerCorrTool_H 1

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "EnergyCorr.h"

/**   
* @class LastLayerCorrTool
*
* Algorithm for calculation of energy leakage correction by correlating
* the energy in the last CAL layer with the total.  
*
*
* $Header$
*/


class LastLayerCorrTool :  public EnergyCorr {

public:
    
    //! destructor
    LastLayerCorrTool( const std::string& type, const std::string& name, const IInterface* parent);
     ~LastLayerCorrTool() {}; 
    
     StatusCode initialize();

        
//! Leakage correction using last layer correlation
/*!This method uses the correlation between the escaping energy and
*  the energy deposited in the last layer of the calorimeter.  Indeed,
* the last layer carries the most important information concerning the leaking
* energy:  the total number of particles escaping through the back should be
* nearly proportional to the  energy deposited in the last layer.
*  The measured signal in that layer can therefore be modified to account
* for the leaking energy. 
* We used the Monte Carlo simulation of the GLAST beam test configuration to 
* determine this correlation at several energies, from 2 GeV up to 40 GeV.  
* For one particular incident energy, the bidimensionnal distribution of 
* the energy escaping and the energy deposited in the last layer can be 
* fitted by a simple linear function:
* \f[ E_{leak} = \alpha E_{last} \f]
* The coefficient \f$\alpha\f$   is parametrized as a function of logarithm
 of incident energy: 
*\f[ \alpha = (p_0 + p_1 * lnE )/(1+exp(-p_2*(lnE - p_3))) \f]
* where coefficients \f$ p_i \f$ are the linear functions
* of \f$ cos(\theta) \f$:
* \f[ p_i = a_i + b_i * cos(\theta) \f]
* The coefficients \f$ a_i \f$ and \f$ b_i \f$ are obtained by fitting MC
* simulation results. 
* The value of \f$ cos(\theta) \f$ is determined from tracker reconstruction
* or, if it is not available, from fitting calorimeter data
* (see Fit_Direction() function for details).  
* The only information on the incident energy we have initially is       
*  the total measured energy
* \f$E_m\f$, thus we have to use it as the estimator of incident energy \f$E_0\f$.
*  The reconstructed energy is then:
* \f[ E_{rec} = E_m + \alpha(E_m) E_{last} \f]
* To improve the result, we make one iteration using the new energy estimator to determine
* the  correct value of \f$\alpha\f$.  
*
*\par The method takes 2 arguments:
*\param eTotal Total energy measured in the calorimeter in MeV
*\param elast  Total energy measured in the last layer in MeV
*
*\return Corrected energy in MeV
*
* \warning Give biased results when the last layer gains are misaligned
*
* \author Regis Terrier
*
* \b Revision: 
* - 10/17/00    RT    comments added
* - 09/15/00    RT    bias correction added, tuned on MC data
* - 08/20/00    RT    first implementation
*/
           
     // worker function to calculate energy correction 
     StatusCode doEnergyCorr(double eTotal, Event::CalCluster* cluster);

     StatusCode finalize();
    
     StatusCode execute();

    
    

};

#endif



