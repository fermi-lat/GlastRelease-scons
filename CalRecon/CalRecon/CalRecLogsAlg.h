

#ifndef __CALRECLOGSALG_H
#define __CALRECLOGSALG_H 1

#include "GaudiKernel/Algorithm.h"
#include "ICalGeometrySvc.h"
class CalPedCalib;
class CalCalibLogs;
class CalCalibLog;

class CalADCLogs;
class CalRecLogs;

class CalRecLog;
class CalADCLog;
class CalDetGeo;

//----------------------------------------------
//
//   CalRecLogsAlg
//
//   Algorithm Data constructor of CalRecLogs
//----------------------------------------------
//             J.A Hernando, Santa Cruz 02/29/00
//----------------------------------------------
//##########################################################
class CalRecLogsAlg : public Algorithm
//##########################################################
{
public:
    
    // constructor
    CalRecLogsAlg(const std::string& name, ISvcLocator* pSvcLocator);
    virtual ~CalRecLogsAlg() {}
    
    // operations
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
protected:
    StatusCode retrieve();
    
private:
    
    void computeEnergy(CalRecLog* recLog, const CalADCLog* adcLog, 
        const CalADCLog* pedLog, const CalCalibLog* calibLog);
    void computePosition(CalRecLog* recLog, const CalDetGeo* geoLog,
        const CalCalibLog* calibLog);
    
private:
    
    ICalGeometrySvc* m_CalGeo;
    CalPedCalib* m_CalPedLogs;
    CalCalibLogs* m_CalCalibLogs;
    CalADCLogs* m_CalRawLogs;
    
    CalRecLogs* m_CalRecLogs;
    
    std::string m_PedFileName;
    std::string m_GainFileName;
    std::string m_IntlinFileName;
    std::string m_RailFileName;
    std::string m_SlopeFileName;
    std::string m_ChargePeaksFileName;
    std::string m_MuPeaksFileName;

};

#endif
