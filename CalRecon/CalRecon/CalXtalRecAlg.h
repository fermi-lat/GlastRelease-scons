#ifndef __CALXTALRECALG_H
#define __CALXTALRECALG_H 1

#include "GaudiKernel/Algorithm.h"
#include "Event/Digi/CalDigi.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"





/** @class CalXtalRecAlg
 *  @brief  Calorimeter crystal reconstruction algorithm
 *
 * This algorithm reconstructs energy and position in each calorimeter crystal
 *
 *  @author           A.Chekhtman
 *
 * $Header$
 */
class CalXtalRecAlg : public Algorithm
{
public:
    
    /// constructor
    CalXtalRecAlg(const std::string& name, ISvcLocator* pSvcLocator);
    virtual ~CalXtalRecAlg() {}
    
    // operations
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
protected:
///  function for setting pointers to the input and output data in Gaudi TDS
    StatusCode retrieve();
    
private:
    /// method to calculate energy deposited in a crystal
    void computeEnergy(Event::CalXtalRecData* recLog,
                       const Event::CalDigi* adcLog);
    
    /// method to calculate longitudinal position in a crystal 
    void computePosition(Event::CalXtalRecData* recLog);
    
private:
    /// pointer to input data collection in TDS
    Event::CalDigiCol* m_CalDigiCol;

    /// pointer to the output data collection in TDS
    Event::CalXtalRecCol* m_CalXtalRecCol;

    /// constants defining the position of the fields in VolumeIdentifier 
    enum {fLATObjects, fTowerY, fTowerX, fTowerObjects, fLayer,
            fMeasure, fCALXtal,fCellCmp, fSegment};
        
        
        // constants defined in xml files
    int m_xNum;  ///< x tower number
    int m_yNum;  ///< y tower number
    int m_nTowers;  ///< total number of towers

    int m_eTowerCAL;
    ///< the value of fTowerObject field, defining calorimeter module 
    int m_eLATTowers; ///< the value of fLATObjects field, defining LAT towers 
    int m_CALnLayer;  ///< number of CAL layers
    int m_nCsIPerLayer;  ///< number of Xtals per layer

    int m_eXtal;      ///< the value of fCellCmp field defining CsI crystal
    int m_nCsISeg;  ///< number of geometric segments per Xtal

    int m_eDiodeMSmall;
    ///< the value of fcellCmp field defining small diode at minus face  

    int m_eDiodePSmall;
    ///< the value of fcellCmp field defining small diode at plus face 

    int m_eDiodeMLarge;
    ///< the value of fcellCmp field defining large diode at minus face 

    int m_eDiodePLarge;
    ///< the value of fcellCmp field defining large diode at plus face 

    int m_eMeasureX;
    ///< the value of fmeasure field for crystals along X direction

    int m_eMeasureY;
    ///< the value of fmeasure field for crystals along Y direction

    int m_ePerMeV[2];  ///< gain - electrons/MeV 0=Small, 1=Large
    int m_noise[2];  ///< noise for diodes 0=Small, 1=Large units=electrons
    int m_pedestal;  ///< single pedestal
    int m_maxAdc;  ///< max value for ADC
    int m_thresh;  ///< zero suppression threshold
    double m_maxEnergy[4];  ///< highest energy for each energy range
    double m_lightAtt;  ///< light attenuation factor
    double m_CsILength;  ///< Xtal length
    IGlastDetSvc* detSvc; ///< pointer to the Glast Detector Service

};


#endif
