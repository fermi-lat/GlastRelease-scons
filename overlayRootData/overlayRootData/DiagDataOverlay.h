#ifndef ROOT_DiagDataOverlay_H
#define ROOT_DiagDataOverlay_H

#include "TObject.h"

#include "TkrDiagDataOverlay.h"
#include "CalDiagDataOverlay.h"
#include <vector>

/*! @class DiagDataOverlay
@brief Storage of Diagnostic EM data (trigger primitives).

Consists of 2 vectors, one for CAl, the other for TKR
*/
class DiagDataOverlay : public TObject
{
public:
    DiagDataOverlay()  
    { 
        Clear(); 
    }
    
    DiagDataOverlay(const DiagDataOverlay& d) 
    {
        m_tkr        = d.m_tkr;
        m_cal        = d.m_cal;
        m_lenInBytes = d.m_lenInBytes;
        m_exist      = d.m_exist;
    }
    ~DiagDataOverlay() {}

    void Clear(Option_t *option="") 
    { 
        m_tkr.clear(); m_cal.clear(); 
        m_lenInBytes = 0; 
        m_exist      = false; 
    }

    void Print(Option_t *option="") const
    {
        if( !exist()) 
        {
            printf("No Diagnostic Data\n");                
            return;
        }
        
        printf("Diagnostic Data:\n");
        printf("Len: %lu\n\n", m_lenInBytes);
        printf("CAL Diagnostics:\n");
        std::vector<CalDiagDataOverlay>::const_iterator it;
        for (it = m_cal.begin(); it!= m_cal.end(); it++) 
        {
            printf("Tower: %d, Layer: %d, DataWord: %u\n",
               it->tower(), it->layer(), it->dataWord());
        }
        
        printf("TKR Diagnostics:\n");
        std::vector<TkrDiagDataOverlay>::const_iterator tkrit;
        for (tkrit = m_tkr.begin(); tkrit!= m_tkr.end(); tkrit++) 
        {
            printf("Tower: %d, GTCC: %d, DataWord: %u\n", 
               tkrit->tower(), tkrit->gtcc(), tkrit->dataWord());
        }
    }

    // Added 03.19.2005 by awb to get all Tkr and Cal diagnostics contributions:
    void addCalDiagnostic(const CalDiagDataOverlay &cal) 
    { 
        m_cal.push_back(cal); 
    }
    
    const CalDiagDataOverlay& getCalDiagnosticByIndex(UInt_t index) const
    { 
        return m_cal[index]; 
    }

    Int_t getNumCalDiagnostic() const 
    { 
        return m_cal.size(); 
    }

    void addTkrDiagnostic(const TkrDiagDataOverlay &tkr) 
    { 
        m_tkr.push_back(tkr); 
    }
    
    const TkrDiagDataOverlay& getTkrDiagnosticByIndex(UInt_t index) const 
    { 
        return m_tkr[index]; 
    }
    
    Int_t getNumTkrDiagnostic() const 
    { 
        return m_tkr.size(); 
    }

   void setExist()      {m_exist = true;}
   Bool_t exist() const {return m_exist;}

   void initLength(ULong_t length) {m_lenInBytes = length;}
   ULong_t lenInBytes() const      {return m_lenInBytes;}
   
   void initPacketError(UInt_t packetError) {m_packetError=packetError;}
   UInt_t packetError() const               {return m_packetError;}

private:

    // Added 03.19.2005 awb to get all Tkr and Cal diagnostics contributions:
    std::vector<CalDiagDataOverlay> m_cal;
    std::vector<TkrDiagDataOverlay> m_tkr;


    Bool_t  m_exist; 

    UInt_t  m_packetError;

    ULong_t m_lenInBytes;

    ClassDef(DiagDataOverlay,1) // Storage for GEM TileList
};

#endif
