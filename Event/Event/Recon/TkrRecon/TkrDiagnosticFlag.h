/** @file TkrDiagnosticFlag.h
*
* @author Bill Atwood, Leon Rochester, Johann Cohen-Tanugi, Tracy Usher
*
* $Header$
*
*/

#ifndef TkrDiagnosticFlag_H
#define TkrDiagnosticFlag_H

#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IInterface.h"

#include "idents/TkrId.h"

static const CLID& CLID_TkrDiagnosticFlag = InterfaceID("TkrDiagnosticFlag",  1, 1);

namespace Event { // Namespace

/** 
* @class TkrDiagnosticFlag
*
* @brief Diagnostics output class
*
*/

class TkrDiagnosticFlag : virtual public DataObject
{
public:
    /// Default (null) constructor (just in case...)
    TkrDiagnosticFlag() : m_diagnosticFlag(1) {};

    /// Construct all but the track parameters, they must be set during recon stage
    TkrDiagnosticFlag(bool diagnosticFlag) :
                       m_diagnosticFlag(diagnosticFlag) {};

    //! Destructor
    virtual ~TkrDiagnosticFlag() {return;}

    //! Gaudi stuff: Retrieve pointer to class defininition structure
    virtual const CLID& clID()                 const   { return TkrDiagnosticFlag::classID(); }
    static  const CLID& classID()                      { return CLID_TkrDiagnosticFlag;       }

    inline const bool getDiagnosticFlag()      const {return m_diagnosticFlag;}
    inline void setDiagnosticFlag(const bool diagnosticFlag) {
        m_diagnosticFlag = diagnosticFlag;
    }

    inline std::ostream& fillStream( std::ostream& s ) const;
    
private:
    /// Define here variables to keep diagnostic information for each event
    bool m_diagnosticFlag;          // 0/1 means don't/do make and use diagnostic info
};

inline std::ostream& Event::TkrDiagnosticFlag::fillStream( std::ostream& s ) const 
{ 
  s << "diagnosticFlag : " << getDiagnosticFlag()         << "\n";

  return s; 
}

}; //Namespace

#endif
