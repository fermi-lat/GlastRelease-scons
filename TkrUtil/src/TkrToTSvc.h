/*
@file TkrToTSvc.h

@brief keeps track of the left-right splits of the tracker planes
@author Leon Rochester

$Header$

*/
#ifndef TkrToTSvc_H
#define TkrToTSvc_H 1

// Include files
#include "TkrUtil/ITkrToTSvc.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "GaudiKernel/Service.h"
#include "idents/TowerId.h"

/** @class TkrToTSvc
* @brief Service to store and compare to a list of desired failure modes in
* the TKR.
*
* Author:  L. Rochester (after R.Dubois)
*
*/

class TkrToTSvc : public Service, virtual public ITkrToTSvc  {

public:

    enum {NTOWERS=16, NLAYERS=18, NVIEWS=2, NSTRIPS=1536};

    TkrToTSvc(const std::string& name, ISvcLocator* pSvcLocator); 

    StatusCode initialize();
    StatusCode finalize();

    /// queryInterface - for implementing a Service this is necessary
    StatusCode queryInterface(const IID& riid, void** ppvUnknown);

    static const InterfaceID& interfaceID() {
        return ITkrToTSvc::interfaceID(); 
    }

    /// return the service type
    const IID& type() const;

    double getGain(int tower, int layer, int view, int strip) const;
    double getGain2(int tower, int layer, int view, int strip) const; 
    double getThreshold(int tower, int layer, int view, int strip) const; 
    double getQuality(int tower, int layer, int view, int strip) const; 


    double getMuonFactor(int tower, int layer, int view, int strip) const 
        {
        if (valid(tower, layer, view, strip)) {
            return (double) m_ToTQuality[tower][layer][view][strip];
        }else {
            return -1.;
        }
    }
    double getCountsPerMicrosecond() const { return m_countsPerMicrosecond;}
    double getMevPerMip() const { return m_mevPerMip; }
    double getFCPerMip() const { return m_fCPerMip; }
    int    getMaxToT() const { return m_maxToT; }

    double getCharge(double ToT, int tower, int layer, int view, int strip) const;
    double getMipsFromToT(double ToT, int tower, int layer, int view, int strip) const;
    double getMipsFromCharge(double charge, int tower, int layer, int view, int strip) const;
    int    getRawToT(double eDep, int tower, int layer, int view, int strip) const;

        /// update the pointer
    void update(CalibData::TkrTotCol* pToT) { m_pToT = pToT; }


private:
    /// internal init method
    StatusCode doInit();

    /// check index
    bool valid(int tower, int layer, int view, int strip) const
    {
        return (tower>-1 && tower <NTOWERS && layer>-1 && layer<NLAYERS
            && view>-1 && view<NVIEWS && strip>-1 && strip<NSTRIPS);
    }

    /// mode: currently "default" or "EM"
    std::string m_mode;
    /// name of file containing splits
    std::string m_ToTFile;
    /// default Gain
    double m_defaultGain;
    /// default quadratic term;
    double m_defaultGain2;
    /// default Threshold
    double m_defaultThreshold;
    /// default quality factor
    double m_defaultQuality;
    /// default muon correction factor
    double m_defaultMuonFactor;
    /// ToT counts per microsecond
    double m_countsPerMicrosecond;
    /// Energy deposited by a mini particle traversing a silicon plane
    double m_mevPerMip;
    /// Charge deposited by a mini particle traversing a silicon plane
    double m_fCPerMip;
    /// array of gains, in microseconds/fC
    int    m_maxToT;
    float m_ToTGain      [NTOWERS][NLAYERS][NVIEWS][NSTRIPS];
    /// array of quadratic terms, in microseconds/fC**2
    float m_ToTGain2     [NTOWERS][NLAYERS][NVIEWS][NSTRIPS];
    /// array of Thresholds, in microseconds = extrapolation to zero charge
    float m_ToTThreshold [NTOWERS][NLAYERS][NVIEWS][NSTRIPS];
    /// array of quality factors
    float m_ToTQuality   [NTOWERS][NLAYERS][NVIEWS][NSTRIPS];
    /// array of Muon correction normalizations, should be order(1)
    float m_ToTMuonFactor [NTOWERS][NLAYERS][NVIEWS][NSTRIPS];
    /// pointer to geometry service
    ITkrGeometrySvc* m_tkrGeom;
    ///
    CalibData::TkrTotCol* m_pToT;
};


#endif // TkrToTSvc_H

