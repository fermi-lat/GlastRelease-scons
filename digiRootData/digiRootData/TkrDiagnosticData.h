#ifndef ROOT_TkrDiagnosticData_H
#define ROOT_TkrDiagnosticData_H

#include "TObject.h"

/** @class TkrDiagnosticData
 * @brief Store the EM TkrDiagnosticData trigger primitives
 *
 * $Header$
*/

class TkrDiagnosticData : public TObject {
public:
            
    TkrDiagnosticData() : m_datum(0) {};
    TkrDiagnosticData(UInt_t datum) : m_datum(datum) {}
    virtual ~TkrDiagnosticData() {}

    UInt_t getDataWord() const { return m_datum; };

    void initialize(UInt_t datum);

    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;

    UInt_t GTRC(Int_t gtrc) const;

private:
    /// packed word containing trigger primitive for these TKR layers
    UInt_t m_datum;    

    ClassDef(TkrDiagnosticData,1) // EM TkrDiagnosticData information
};

#endif
