#ifndef LDFCONVERTER_DIAGNOSTICDATA_H
#define LDFCONVERTER_DIAGNOSTICDATA_H

#include <iostream>
#include <vector>

#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/StreamBuffer.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/ObjectVector.h"

#include "Event/TopLevel/Definitions.h"
#include "Event/TopLevel/EventModel.h"

// From EBF Online library
#include "CALdiagnostic.h"
#include "TKRdiagnostic.h"

/**
* @class DiagnosticData
* @brief TDS for storing the diagnostic (trigger primitives)

* Constains 2 local classes for CAL and TKR trigger primitives
*/
extern const CLID& CLID_LdfDiagnosticData;
namespace LdfEvent {
    class CalDiagnosticData {
    public:
        CalDiagnosticData(unsigned datum) 
            : m_datum(datum)
        {
        }

        CalDiagnosticData()  {}
        ~CalDiagnosticData() {}

        unsigned dataWord() const { return m_datum; };

    private:
        unsigned m_datum;
    };

    class TkrDiagnosticData {
    public:

        TkrDiagnosticData(unsigned datum) 
            : m_datum(datum) {}
            TkrDiagnosticData()  {}
            ~TkrDiagnosticData() {}

            unsigned dataWord() const { return m_datum; };

    private:
        unsigned m_datum;
    };

    class DiagnosticData : public DataObject{
    public:
        virtual ~DiagnosticData();
        virtual std::ostream& fillStream(std::ostream &s) const;
        friend std::ostream& operator << (std::ostream &s, const DiagnosticData &obj);

        DiagnosticData();

        inline void addCalDiagnostic(const CalDiagnosticData& cal) { m_cal.push_back(cal); };
        inline int getNumCalDiagnostic() const { return m_cal.size(); };
        const CalDiagnosticData& getCalDiagnosticByIndex(unsigned int ind)const { return m_cal[ind]; };
        inline void addTkrDiagnostic(const TkrDiagnosticData& tkr) { m_tkr.push_back(tkr); };
        inline int getNumTkrDiagnostic() const { return m_tkr.size(); };
        const TkrDiagnosticData& getTkrDiagnosticByIndex(unsigned int ind)const { return m_tkr[ind]; };

    private:
        std::vector<CalDiagnosticData> m_cal;
        std::vector<TkrDiagnosticData> m_tkr;
    };

    // inline the public get methods for clients besides OnboardFilter.

    // the put methods, here for now
    inline DiagnosticData::DiagnosticData(){
        m_cal.clear();
        m_tkr.clear();
    }

    inline DiagnosticData::~DiagnosticData(){
    }


    inline std::ostream& DiagnosticData::fillStream(std::ostream &s) const{
        int i;
        for (i = 0; i < getNumCalDiagnostic(); i++) {
            s<<"Full CAL Data Word # " << i << " " << m_cal[i].dataWord() <<std::endl;
        }
        for (i = 0; i < getNumTkrDiagnostic(); i++) {
            s<<"Full TKR Data Word # " << i << " " <<m_tkr[i].dataWord()<<std::endl;
        }
        return s;
    }


    inline std::ostream& operator<<(std::ostream &s, const DiagnosticData &obj){
        return obj.fillStream(s);
    }

}//namespace LdfEvent 


#endif
