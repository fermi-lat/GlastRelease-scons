#ifndef LDFCONVERTER_ERRORDATA_H
#define LDFCONVERTER_ERRORDATA_H

#include <iostream>
#include <vector>

#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/StreamBuffer.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/ObjectVector.h"

#include "Event/TopLevel/Definitions.h"
#include "Event/TopLevel/EventModel.h"

/**
* @class ErrorData
* @brief TDS for storing the error
*/

static const CLID& CLID_LdfErrorData = InterfaceID("LdfErrorData", 1, 0);

namespace LdfEvent {
    class TowerErrorData {
    public:
        TowerErrorData();
        TowerErrorData(unsigned short tower, unsigned short cal, unsigned short tkr, bool phs, bool tmo);

        virtual ~TowerErrorData();
        void clear();
        void init(unsigned short tower, unsigned short cal, unsigned short tkr, bool phs, bool tmo);

        unsigned short tower() const { return m_tower; }
        unsigned short cal() const { return m_cal; }
        unsigned short tkr() const { return m_tkr; }
        bool phs() const { return m_phs; }
        bool tmo() const { return m_tmo; }

    private:
        unsigned short m_tower;
        unsigned short m_cal;  // CAL
        unsigned short m_tkr;  // TKR
        bool m_phs;  // Phasing Error
        bool m_tmo;  // Timeout Error
    };

    class ErrorData : public DataObject {
    public:
        ErrorData() { clear(); };
        virtual ~ErrorData() { clear(); };
        void clear() { m_errorCol.clear(); };


       inline void addTowerError(const TowerErrorData& err) { 
                                     m_errorCol.push_back(err); };
        inline int getNumTowerError() const { return m_errorCol.size(); };
        const TowerErrorData& getTowerErrorByIndex(unsigned int ind)const
{ return m_errorCol[ind]; };

        const std::vector<TowerErrorData>& errorCol() const { 
                                                      return m_errorCol; };

        virtual std::ostream& fillStream(std::ostream &s) const;
        friend std::ostream& operator << (std::ostream &s, const ErrorData &obj);
        /// Retrieve reference to class definition structure
        virtual const CLID& clID() const   { return ErrorData::classID(); }
        static const CLID& classID()       { return CLID_LdfErrorData; };

    private:
        std::vector<TowerErrorData> m_errorCol;
    };


    // the put methods, here for now
    inline TowerErrorData::TowerErrorData(){
        clear();
    }

    inline TowerErrorData::TowerErrorData(unsigned short tower, 
                                        unsigned short cal, unsigned short tkr,
                                        bool phs, bool tmo) {
        init(tower, cal, tkr, phs, tmo);
    }

    inline void TowerErrorData::init(unsigned short tower, unsigned short cal, unsigned short tkr, bool phs, bool tmo) {
        clear();
        m_cal = cal;
        m_tkr = tkr;
        m_phs = phs;
        m_tmo = tmo;
        m_tower = tower;
    }

    inline TowerErrorData::~TowerErrorData(){
    }

    inline void TowerErrorData::clear() {
        m_tower = 0;
        m_cal = 0;
        m_tkr = 0;
        m_phs = false;
        m_tmo = false;
    }





    inline std::ostream& ErrorData::fillStream(std::ostream &s) const{
        int i;
        for (i = 0; i < getNumTowerError(); i++) {
            s<<"Error " << "CAL: " << m_errorCol[i].cal() << " TKR: " 
             << m_errorCol[i].tkr() << " Phase: " << m_errorCol[i].phs()
             << " TimeOut: " << m_errorCol[i].tmo() <<std::endl;
        }
        return s;
    }


    inline std::ostream& operator<<(std::ostream &s, const ErrorData &obj){
        return obj.fillStream(s);
    }

}//namespace LdfEvent 


#endif
