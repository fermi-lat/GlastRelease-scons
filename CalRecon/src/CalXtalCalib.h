#ifndef CALXTALCALIB_H
#define CALXTALCALIB_H

#include <iostream>
#include <map>
#include "idents/CalXtalId.h"


template <class CalibData, class CalibElement, unsigned int nFaces=2>
class CalXtalCalib
{
public:

    CalXtalCalib(){}

    ~CalXtalCalib(){}

    const CalibData* getRangeCalib( unsigned short range,
                             idents::CalXtalId::XtalFace face
                                  =idents::CalXtalId::POS) const;
    void addRangeCalib(CalibData* rangeCalib, unsigned int range,
                         idents::CalXtalId::XtalFace face = idents::CalXtalId::POS);

    void addElement(CalibElement* elem, const unsigned range,
                        const idents::CalXtalId::XtalFace face
                             =idents::CalXtalId::POS);

    void generateCalib(unsigned int nRanges);
        
private:

    std::map<int,CalibData*> m_calib[nFaces];

};


template <class CalibData, class CalibElement, unsigned int nFaces>
const CalibData* CalXtalCalib<CalibData,CalibElement,nFaces>::getRangeCalib(unsigned short range,
                              idents::CalXtalId::XtalFace face
                                  =idents::CalXtalId::POS) const
    {
        const std::map<int,CalibData*>& calibFace = 
        (nFaces>1 && face == idents::CalXtalId::NEG) ? m_calib[1]:m_calib[0];

        return (*(calibFace.find(range))).second; 

    }

template <class CalibData, class CalibElement, unsigned int nFaces>
    void CalXtalCalib<CalibData,CalibElement,nFaces>::addRangeCalib(CalibData* rangeCalib, unsigned int range,
                        idents::CalXtalId::XtalFace face)


    {
        
        std::map<int,CalibData*>& calibFace = 
        (nFaces>1 && face == idents::CalXtalId::NEG) ? m_calib[1]:m_calib[0];

        calibFace[range] = rangeCalib;

    }

template <class CalibData, class CalibElement, unsigned int nFaces>
    void CalXtalCalib<CalibData,CalibElement,nFaces>::addElement(CalibElement* elem, const unsigned range,
                        const idents::CalXtalId::XtalFace face
                             =idents::CalXtalId::POS) 
    {
        std::map<int,CalibData*>& calibFace = 
        (nFaces>1 && face == idents::CalXtalId::NEG) ? m_calib[1]:m_calib[0];
           

        (calibFace[range])->addElement(elem);
    }

template <class CalibData, class CalibElement, unsigned int nFaces>
        void CalXtalCalib<CalibData,CalibElement,nFaces>::generateCalib(unsigned int nRanges)
    {


        unsigned int  range, iface;

                    for (range=0; range < nRanges; range++)
                    {
                        for (iface=0; iface < nFaces; iface++)
                        {

                            idents::CalXtalId::XtalFace face = 
                                iface==0 ? idents::CalXtalId::POS : idents::CalXtalId::NEG;
                            
                            CalibData* rangeCalib = new CalibData();
            
                            rangeCalib->generateCalib();

                            addRangeCalib(rangeCalib,range,face);                            
                        }
                    }

    }


#endif // CALXTALCALIB_H 
