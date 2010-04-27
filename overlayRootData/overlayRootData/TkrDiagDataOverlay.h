#ifndef ROOT_TkrDiagDataOverlay_H
#define ROOT_TkrDiagDataOverlay_H

#include "TObject.h"

/* @class TkrDiagDataOverlay
@brief Stores TKR trigger primitives
*/

class TkrDiagDataOverlay : public TObject
{
public:

    TkrDiagDataOverlay(UInt_t datum, UInt_t tower, UInt_t gtcc) : m_datum(datum),m_tower(tower),m_gtcc(gtcc) {}

    TkrDiagDataOverlay(const TkrDiagDataOverlay &t) 
    {
        m_datum = t.m_datum;
        m_tower = t.m_tower;
        m_gtcc = t.m_gtcc;
    }

    TkrDiagDataOverlay()  
    {
        Clear();
    }
    
    ~TkrDiagDataOverlay() {}

    void Clear(Option_t *option="") 
    {
        m_datum = 0;
        m_tower = 0;
        m_gtcc = 0;
    }

    void initialize(UInt_t data, UInt_t tower, UInt_t gtcc) 
    { 
         m_datum = data; 
         m_tower=tower; 
         m_gtcc=gtcc; 
    }

    UInt_t dataWord() const { return m_datum;}
    UInt_t tower()    const { return m_tower;}
    UInt_t gtcc()     const { return m_gtcc;}

    void Print(Option_t *option="") const {}

private:
    UInt_t m_datum;
    UInt_t m_tower;
    UInt_t m_gtcc;

    ClassDef(TkrDiagDataOverlay,1) // Storage for GEM TileList
};

#endif
