#ifndef ROOT_CalDiagnosticData_H
#define ROOT_CalDiagnosticData_H

#include "TObject.h"

/** @class CalDiagnosticData
 * @brief Store the EM CalDiagnosticData trigger primitives
 *
 * $Header$
*/

class CalDiagnosticData : public TObject {
public:
            
    CalDiagnosticData();
    CalDiagnosticData(UInt_t datum);
    virtual ~CalDiagnosticData();

    void initialize(UInt_t datum);

    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;

    UInt_t getDataWord() const { return m_datum; };
    UInt_t high(UInt_t sign) const;
    UInt_t low(UInt_t sign) const;
    UInt_t logAccepts(UInt_t sign) const;


private:
    /// packed word containing trigger primitives for this set of logs
    UInt_t m_datum;

    ClassDef(CalDiagnosticData,1) // EM CalDiagnosticData information
};

#endif
