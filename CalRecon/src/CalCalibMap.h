#ifndef CALCALIBMAP_H
#define CALCALIBMAP_H

#include <fstream>
#include <map>
#include "idents/CalXtalId.h"
#include "CalXtalCalib.h"

template <class CalibData, class CalibElement, unsigned int nFaces=2>
class CalCalibMap
{
public:

    typedef CalXtalCalib<CalibData,CalibElement,nFaces> XTALCALIB;

    CalCalibMap() {};
       
    ~CalCalibMap() {};

    const XTALCALIB* 
        getXtalCalib(const idents::CalXtalId xtalId) const;


    // initializing calibration data for the crystal with given xtalId
    void addElement(CalibElement* elem, const idents::CalXtalId xtalId,
        unsigned int range, idents::CalXtalId::XtalFace face =
        idents::CalXtalId::POS);

    // reading calibration data from ascii file

       void readFile(std::ifstream& file);


        void writeFile(std::ofstream &file);
        
        void addRangeCalib(CalibData* rangeCalib, idents::CalXtalId xtalId,
            unsigned int range, idents::CalXtalId::XtalFace face);

        void generateCalib(unsigned int nTowers,
                                        unsigned int nLayers,
                                        unsigned int nColumns,
                                        unsigned int nRanges);
        void addXtalCalib(XTALCALIB* xtalCalib, idents::CalXtalId xtalId);
private:

    std::map<idents::CalXtalId,XTALCALIB*> m_calibMap;

};



template <class CalibData, class CalibElement,unsigned int nFaces>
const CalCalibMap<CalibData,CalibElement,nFaces>::XTALCALIB* 
        CalCalibMap<CalibData,CalibElement,nFaces>::getXtalCalib(const idents::CalXtalId xtalId) const
        { return (*(m_calibMap.find(xtalId))).second;}


    // initializing calibration data for the crystal with given xtalId
template <class CalibData, class CalibElement,unsigned int nFaces>
    void CalCalibMap<CalibData,CalibElement,nFaces>::addElement(CalibElement* elem, const idents::CalXtalId xtalId,
        unsigned int range, idents::CalXtalId::XtalFace face =
        idents::CalXtalId::POS)
    {(m_calibMap[xtalId])->addElement(elem,range,face);}

    // reading calibration data from ascii file

template <class CalibData, class CalibElement, unsigned int nFaces>
       void CalCalibMap<CalibData, CalibElement,nFaces>::readFile(std::ifstream& file)
    {

        unsigned short tower, layer, col, range, iface;
        
        while ((file >> tower >> layer >> col >> range).good())
        {
            const idents::CalXtalId xtalId(tower,layer,column);

            iface=0; if(nFaces>1) file >> iface;
            idents::CalXtalId::XtalFace face = 
                iface==0 ? idents::CalXtalId::POS : idents::CalXtaId::NEG;

            CalibData::Element* elem = new CalibData::Element();
            
            elem->readFile(file);

            addElement(elem,xtalId,range,face);
        }


    }


template <class CalibData, class CalibElement,unsigned int nFaces>
        void CalCalibMap<CalibData, CalibElement,nFaces>::writeFile(std::ofstream &file)
    {
        unsigned int tower,layer,col,range,iface;
        std::map<idents::CalXtalId,XTALCALIB*>::iterator itXtal;

        for(itXtal=m_calibMap.begin();itXtal!==m_calibMap.end();itXtal++)
        {
            idents::CalXtalId xtalId = (*itXtal).first();
            tower = xtalId.getTower();
            layer = xtalId.getLayer();
            col   = xtalId.getColumn();

            for(iface = 0; iface < nFaces;iface++)
            {
                idents::CalXtalId::XtalFace face = 
                iface==0 ? idents::CalXtalId::POS : idents::CalXtaId::NEG;
                for (range = 0; range<4; range++)
                {
                    const CalibData* rangeCalib 
                        = (*itXtal).second()->getRangeCalib(range,face);

                    rangeCalib->writeFile(file,tower,layer,col,range,face);
                }
            }
        }

    }

        
template <class CalibData, class CalibElement,unsigned int nFaces>
        void CalCalibMap<CalibData,CalibElement,nFaces>::addRangeCalib(CalibData* rangeCalib, idents::CalXtalId xtalId,
            unsigned int range, idents::CalXtalId::XtalFace face)

        {

            XTALCALIB* xtalCalib = m_calibMap[xtalId];
            xtalCalib->addRangeCalib(rangeCalib,range,face);


        }

template <class CalibData, class CalibElement,unsigned int nFaces>
        void CalCalibMap<CalibData,CalibElement,nFaces>::generateCalib(unsigned int nTowers,
                                        unsigned int nLayers,
                                        unsigned int nColumns,
                                        unsigned int nRanges)
    {


        unsigned int tower, layer, col, range, iface;
        for (tower=0; tower < nTowers; tower++)
        {
            for (layer=0; layer < nLayers; layer++)
            {
                for (col=0; col < nColumns; col++)
                {
                    idents::CalXtalId xtalId(tower,layer,col);

                    addXtalCalib(new XTALCALIB, xtalId);

                    for (range=0; range < nRanges; range++)
                    {
                        for (iface=0; iface < nFaces; iface++)
                        {

                            idents::CalXtalId::XtalFace face = 
                                iface==0 ? idents::CalXtalId::POS : idents::CalXtalId::NEG;
                            
                            CalibData* rangeCalib = new CalibData();
            
                            rangeCalib->generateCalib();
                            
                            
                            addRangeCalib(rangeCalib,xtalId,range,face);


                        }
                    }
                }
            }
        }
    }

template <class CalibData, class CalibElement, unsigned int nFaces>
        void CalCalibMap<CalibData,CalibElement,nFaces>::addXtalCalib(XTALCALIB* xtalCalib, idents::CalXtalId xtalId)
        {

            m_calibMap[xtalId] = xtalCalib;

        }
    


#endif // CALCALIBMAP_H 
