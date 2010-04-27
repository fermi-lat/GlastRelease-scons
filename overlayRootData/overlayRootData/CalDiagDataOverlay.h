#ifndef ROOT_CalDiagDataOverlay_H
#define ROOT_CalDiagDataOverlay_H

#include "TObject.h"

/*
*  @class CalDiagDataOverlay
*  @brief Stores trigger primitives for CAL
*/
class CalDiagDataOverlay : public TObject
{
public:
    CalDiagDataOverlay(UInt_t datum, UInt_t tower, UInt_t layer) 
        : m_datum(datum),m_tower(tower),m_layer(layer) {}

    CalDiagDataOverlay()  {Clear();}
    CalDiagDataOverlay(const CalDiagDataOverlay &c) 
    { 
        m_datum = c.m_datum;
        m_tower = c.m_tower;
        m_layer = c.m_layer;
    }

    ~CalDiagDataOverlay() {}

    void Clear(Option_t *option="") 
    {
        m_datum = 0;
        m_tower = 0;
        m_layer = 0;  
    }

    void initialize(UInt_t data, UInt_t tower, UInt_t layer) 
    { 
         m_datum = data; 
         m_tower = tower; 
         m_layer = layer; 
    }

    UInt_t dataWord() const { return m_datum; }
    UInt_t tower() const { return m_tower; }
    UInt_t layer() const { return m_layer; }

    void Print(Option_t *option="") const {}

private:
    UInt_t m_datum;
    UInt_t m_tower;
    UInt_t m_layer;

    ClassDef(CalDiagDataOverlay,1) // Storage for CAL diagnostic data
};

#endif
