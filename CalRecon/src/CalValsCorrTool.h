
/** @file CalValsCorrTool.h
@brief declaration of the class

$Header$

*/
#ifndef __CalValsCorrTool_H
#define __CalValsCorrTool_H 1

#include "EnergyCorr.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

class IPropagatorSvc;
class ITkrGeometrySvc;
class IPropagator;
class IDataProviderSvc;

/**   
* @class CalValsCorrTool
* @author Bill Atwood
*
* Tool to corrected energy for cracks and leakage.
*
* Copied by THB from AnalysisNtuple::CalValsTool.cxx revision 1.43
*
* $Header$
*/


class CalValsCorrTool : public EnergyCorr {

public:

    //! destructor
    CalValsCorrTool( const std::string& type, const std::string& name, const IInterface* parent);
    ~CalValsCorrTool() {}; 

    StatusCode initialize();

    // worker function to get the corrected energy      
    StatusCode doEnergyCorr( const CalClusteringData *, Event::CalCluster * ) ;

private:

    /// Bill's calculation here
    StatusCode calculate( const CalClusteringData * );
    double activeDist(Point pos, int& view) const;
    double containedFraction(Point pos, double gap, 
        double r, double costh, double phi) const;
    StatusCode aveRadLens(const CalClusteringData * data, Point x0, Vector t0, double radius, int numSamples);

    /// TkrGeometrySvc used for access to tracker geometry info
    ITkrGeometrySvc* m_tkrGeom;

    /// some Geometry
    double m_towerPitch;
    int    m_xNum;
    int    m_yNum;
    
    /// gets the CAL info from detModel
    StatusCode getCalInfo();

    /// CAL vars
    double m_calXWidth;
    double m_calYWidth;
    double m_calZTop;
    double m_calZBot;

    // Internal Variables
    double m_radLen_CsI, m_rms_RL_CsI;
    double m_radLen_Stuff, m_rms_RL_Stuff;
    double m_radLen_Cntr, m_rms_RL_Cntr; 
    double m_radLen_CntrStuff, m_rms_RL_CntrStuff;

    double m_arcLen_CsI; 
    double m_arcLen_Stuff; 
    double m_arcLen_Cntr;  

    IPropagatorSvc* m_propSvc;
    IPropagator * m_G4PropTool; 



    //Global Calorimeter Tuple Items -- but only CAL_Energy_Corr is actually used
    double CAL_EnergySum; 
    double CAL_Leak_Corr;
    double CAL_Edge_Corr; 
    double CAL_EdgeSum_Corr;     
    //double CAL_Total_Corr; 
    double CAL_TotSum_Corr; 
    double CAL_Energy_LLCorr; 

    double CAL_CsI_RLn;
    double CAL_Tot_RLn;
    double CAL_Cnt_RLn; 
    double CAL_LAT_RLn; 
    double CAL_DeadTot_Rat;
    double CAL_DeadCnt_Rat; 
    //double CAL_a_Parm;
    //double CAL_b_Parm; 
    double CAL_t_Pred; 
    double CAL_deltaT;

    double CAL_EneSum_Corr;
    double CAL_Energy_Corr;
    double CAL_xEcntr;
    double CAL_yEcntr;
    double CAL_zEcntr;
    double CAL_xdir;
    double CAL_ydir;
    double CAL_zdir;

    double CAL_TwrEdgeCntr;
    double CAL_TwrEdge0;
    double CAL_LATEdge; 

    double CAL_TE_Nrm;
    double CAL_Track_Sep;

    double CAL_Lyr0_Ratio;
    double CAL_Lyr7_Ratio;
    double CAL_BkHalf_Ratio;

    double CAL_Xtal_Ratio;
    double CAL_Xtal_maxEne; 
    double CAL_eLayer[8];
    double CAL_No_Xtals_Trunc;
    double CAL_Long_Rms;
    double CAL_Trans_Rms;
    double CAL_LRms_Ratio;

    double CAL_MIP_Diff; 
    double CAL_MIP_Ratio;

    //Calimeter items with Recon - Tracks
    double CAL_Track_DOCA;
    double CAL_Track_Angle;
    double CAL_TwrGap; 
    double CAL_x0;
    double CAL_y0;
    double CAL_z0;


};

#endif



