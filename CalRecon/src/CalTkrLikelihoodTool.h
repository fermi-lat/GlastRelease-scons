#ifndef __CalTkrLikelihoodTool_H
#define __CalTkrLikelihoodTool_H 1

#include "LikelihoodTool.h"

/**   
* @class CalTkrLikelihoodTool
*
* Algorithm for correction of energy degradation in the tracker by correlating
* the energy in the CAL with the number of hit TKR strips.  
*
*
*/


class CalTkrLikelihoodTool :  public LikelihoodTool {
  public:
    //! destructor
    CalTkrLikelihoodTool( const std::string& type, const std::string& name, 
                 const IInterface* parent);
    ~CalTkrLikelihoodTool(){}
    
    StatusCode initialize();

  //! Tracker energy degradation correction using number of TKR hit strips
  /*! This method uses the correlation between the energy \"lost\" in the 
  * tracker and the energy deposited in the calorimeter.
  * We used the Monte Carlo simulation of the LAT to determine this correlation
  * at several energies, from 50 MeV up to 1 GeV, and angles from 0 to 32\deg. 
  * For one particular incident energy and angle, the bidimensionnal
  * distribution of  the  number of hit strips and the energy deposited in the 
  * CAL can be characterised by the 1D distribution:
  * \f[ E_{CAL} + \alpha TkrHits \f]
  * where  \f$\alpha\f$ is been optimised so as to obtain the narrowest such
  * distribution, normalised to a probability and with its MPV at the incident
  * energy.
  * These distributions can be used to defined a probability density function.
  * The reconstructed energy for a given event then becomes the one maximising
  * the probability, for a reconstruced direction, CAL raw energy,...
  *
  * \par The method takes 4 arguments:
  * \param eTotal Total energy measured in the calorimeter in MeV
  * \param nHits  Total number of hit strips in the CAL
  * \param vertex[3] reconstructed vertex position  
  * \param dir[3] reconstructed direction
  *
  *\return Corrected energy in MeV
  *
  *\warning needs TKR reconstruction
  *
  *\author
  */
           
     StatusCode doEnergyCorr( const CalClusteringData* , Event::CalCluster* );

    
 private:
    double integratedDistance2TowerSide(double, double) const;

    // CAL origin, tower pitch, and gap between a CAL and its tower,
    // needed for the cuts
    double m_zOriginCAL;
    double m_pitchTOWER;
    double m_halfwidthCAL;
    double m_heightCAL;
    
    /// input data file containing log normal parameters
    std::string	m_dataFile;
};
#endif
