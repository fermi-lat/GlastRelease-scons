/** @file CalClustersAlg.h
    @brief Defintion of CalClustersAlg

*/

#ifndef __CALCLUSTERSALG_H
#define __CALCLUSTERSALG_H 1

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/Property.h"
#include "CalIClusteringTool.h"
#include "IEnergyCorr.h"

class IGlastDetSvc;

/**   
* @class CalClustersAlg
*
* @brief Algorithm for reconstruction of energy and direction of incident particle
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

        
/*!Performs the reconstruction, creates one CalCluster object and stores
 * there the following results: 
 * - Energy per layer is computed and stored in CalCluster in MeV
 * - Barycenter per layer is also computed and stored in CalCluster
 * - The total energy, position and direction for whole calorimeter
 * - Performs high energy corrections: see Profile() and Leak() for details
 */        
    StatusCode execute();

/*! Finalization of algorithm 
 *  - minuit deleted
 */    
    StatusCode finalize();
        
protected:
    
    StatusCode retrieve();
    
private:
    
    //! reconstructed data for crystals, the input of the reconstruction
    Event::CalXtalRecCol* m_calXtalRecCol;
    //! the clusters list, the output of the reconstruction
    Event::CalClusterCol* m_calClusterCol;
    //! the minimizer for Profile()

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

    /// name of Tool for finding clusters
    StringProperty m_clusteringToolName;

    /// pointer to actual tool for finding clusters
    CalIClusteringTool* m_clusteringTool;

    /// name of Tool for finding last layer energy leakage
    StringProperty m_lastLayerToolName;

    /// pointer to actual tool for last layer energy correlation
    IEnergyCorr* m_lastLayerTool;

    /// name of Tool for finding last layer energy leakage
    StringProperty m_profileToolName;

    /// pointer to actual tool for last layer energy correlation
    IEnergyCorr* m_profileTool;

    /// name of Tool for calling CalVals tool
    StringProperty m_calValsCorrToolName;

    /// pointer to actual tool for calling CalVals tool
    IEnergyCorr* m_calValsCorrTool;
};

#endif



