#ifndef __CALXTALRECALG_H
#define __CALXTALRECALG_H 1

#include "GaudiKernel/Algorithm.h"
#include "GlastEvent/Digi/CalDigi.h"
#include "CalRecon/CalXtalRecData.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

class CalDigi;
class CalXtalRecData;

//----------------------------------------------
//
//   CalXtalRecAlg
//
//   Algorithm reconstructing energy and position for each calorimeter crystal
//   ideal case: no noise, non-linearity correction etc.
//
//----------------------------------------------
//             A.Chekhtman, Apr,18, 2002
//----------------------------------------------
//##########################################################
class CalXtalRecAlg : public Algorithm
//##########################################################
{
public:
    
    // constructor
    CalXtalRecAlg(const std::string& name, ISvcLocator* pSvcLocator);
    virtual ~CalXtalRecAlg() {}
    
    // operations
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
protected:
    StatusCode retrieve();
    
private:
    
    void computeEnergy(CalXtalRecData* recLog, const cal::CalDigi* adcLog); 
    void computePosition(CalXtalRecData* recLog);
    
private:
    
//    ICalGeometrySvc* m_CalGeo;
	cal::CalDigiCol* m_CalDigiCol;
    CalXtalRecCol* m_CalXtalRecCol;

        enum {fLATObjects, fTowerY, fTowerX, fTowerObjects, fLayer,
            fMeasure, fCALXtal,fCellCmp, fSegment};
        
        
        /// constants defined in xml files
        int m_xNum;  // x tower number
        int m_yNum;  // y tower number
        int m_nTowers;  // total number of towers
        int m_eTowerCAL;  
        int m_eLATTowers;
        int m_CALnLayer;  // number of CAL layers
        int m_nCsIPerLayer;  // number of Xtals per layer
        int m_eXtal;
        int m_nCsISeg;  // number of geometric segments per Xtal
        int m_eDiodeMSmall;   
        int m_eDiodePSmall;
        int m_eDiodeMLarge;
        int m_eDiodePLarge;
        int m_eMeasureX;
        int m_eMeasureY;
        int m_ePerMeV[2];  // gain - electrons/MeV 0=Small, 1=Large
        int m_noise[2];  // noise for diodes 0=Small, 1=Large units=electrons
        int m_pedestal;  // single pedestal
        int m_maxAdc;  // max value for ADC
        int m_thresh;  // zero suppression threshold
        double m_maxEnergy[4];  // highest energy for each energy range
        double m_lightAtt;  // light attenuation factor
        double m_CsILength;  // Xtal length
	    IGlastDetSvc* detSvc;

};

#endif
