#ifndef ROOT_TkrDiagnosticData_H
#define ROOT_TkrDiagnosticData_H

#include "TObject.h"

/** @class TkrDiagnosticData
 * @brief Store the EM TkrDiagnosticData trigger primitives
 *
 * For details concering the DIAGNOSTIC contribution in the LDF please refer to:
 * http://www-glast.slac.stanford.edu/IntegrationTest/ONLINE/docs/TEM.pdf
 * $Header$
*/

class TkrDiagnosticData : public TObject {
public:
            
    TkrDiagnosticData() : m_datum(0) {};
    TkrDiagnosticData(UInt_t datum) : m_datum(datum) {}
    virtual ~TkrDiagnosticData() {}

    void initialize(UInt_t datum);

    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;

	/**@defgroup TkrDiagnosticDataGroup TkrDiagnosticData End-User Interface */
	/*@{*/
	/// Returns the full packed data word for the TkrDiagnostic
    UInt_t getDataWord() const { return m_datum; };
    UInt_t GTRC(Int_t gtrc) const;
	/*@}*/

private:
    /// packed word containing trigger primitive for these TKR layers
    UInt_t m_datum;    

    ClassDef(TkrDiagnosticData,1) // EM TkrDiagnosticData information
};

#endif
