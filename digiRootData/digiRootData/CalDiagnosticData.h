#ifndef ROOT_CalDiagnosticData_H
#define ROOT_CalDiagnosticData_H

#include "TObject.h"

/** @class CalDiagnosticData
 * @brief Store the EM CalDiagnosticData trigger primitives
 *
 * For details concerning the DIAGNOSTIC contribution in the LDF please refer to:
 * http://www-glast.slac.stanford.edu/IntegrationTest/ONLINE/docs/TEM.pdf
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

	/** @defgroup CalDiagnosticDataGroup CalDiagnosticData End-User Interface */
    /*@{*/
	/// Returns the full packed data word for this CalDiagnosticData object
    UInt_t getDataWord() const { return m_datum; };
    UInt_t high(UInt_t sign) const;
    UInt_t low(UInt_t sign) const;
    UInt_t logAccepts(UInt_t sign) const;
	/*@}*/

private:
    /// packed word containing trigger primitives for this set of logs
    UInt_t m_datum;

    ClassDef(CalDiagnosticData,1) // EM CalDiagnosticData information
};

#endif
