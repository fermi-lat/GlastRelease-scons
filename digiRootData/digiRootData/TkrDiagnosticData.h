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
            
    TkrDiagnosticData() : m_datum(0), m_tower(0), m_gtcc(0) {};
    TkrDiagnosticData(UInt_t datum,UShort_t tower, UShort_t gtcc) 
                     : m_datum(datum), m_tower(tower), m_gtcc(gtcc) {}
    virtual ~TkrDiagnosticData() {}

    void initialize(UInt_t datum, UShort_t tower, UShort_t gtcc);

    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;

	/**@defgroup TkrDiagnosticDataGroup TkrDiagnosticData End-User Interface */
	/*@{*/
	/// Returns the full packed data word for the TkrDiagnostic
    UInt_t getDataWord() const { return m_datum; };
    UInt_t GTRC(Int_t gtrc) const;
    UShort_t tower() const {return m_tower; };
    UShort_t gtcc() const { return m_gtcc; };
	/*@}*/

private:
    /// packed word containing trigger primitive for these TKR layers
    UInt_t m_datum;    
    UShort_t m_tower; 
    UShort_t m_gtcc;

    ClassDef(TkrDiagnosticData,2) // EM TkrDiagnosticData information
};

#endif
