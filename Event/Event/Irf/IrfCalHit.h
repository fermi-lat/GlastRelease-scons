#ifndef IrfCalHit_H
#define IrfCalHit_H 1


// Include files
#include <iostream>
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/SmartRefVector.h"
#include "GlastEvent/TopLevel/Definitions.h"
#include "GlastEvent/Utilities/CellID.h"
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/ObjectList.h"


//! Externals 
extern const CLID& CLID_IrfCalHit;
/*! Represents the data from the CsI logs. Originally adapted from H. Arrighi
    IrfAcdHit class. The IRF parser used with this class is CsIDetector
*/



//------------------------------------------------------------------------------
//
// GlastName:   IrfCalHits
//
// Notes:       This class was reworked from the class of the same name in  
//              LHCbEvent
//
//
//------------------------------------------------------------------------------

class IrfCalHit : virtual public ContainedObject {

public:
  //! Constructors
  IrfCalHit()                                    { }
  //! Destructor
  virtual ~IrfCalHit()                           { }

  //! Retrieve pointer to class defininition structure
  virtual const CLID& clID() const   { return IrfCalHit::classID(); }
  static const CLID& classID()       { return CLID_IrfCalHit; }

  //! Retrieve energy
  double energy() const                                 { return m_energy; }
  //! Update energy
  void setEnergy( double value )                        { m_energy = value; }


  //! Serialize the object for reading
  virtual StreamBuffer& serialize( StreamBuffer& s );
  //! Serialize the object for writing
  virtual StreamBuffer& serialize( StreamBuffer& s ) const;
  //! Fill the ASCII output stream
  virtual std::ostream& fillStream( std::ostream& s ) const;

  
  //! response from the plus end of the log
  double plusResponse() const               { return m_plusEnd; }
  void  setPlusResponse(float value)       { m_plusEnd = value; }

  //! reponse from the minus end of the log
  double minusResponse() const                { return m_minusEnd; }
  void  setMinusResponse(float value)        {m_minusEnd = value;}

  //!  The Layer
  unsigned int layer () const               { return m_layer; }
  void setLayer(unsigned int value)         { m_layer = value; }

  //! Log index number
  unsigned int index () const               { return m_index; }
  void setId (long value)                   { m_id = value; }



private:
  
  double m_energy;  //! Energy
  double m_minusEnd;  // response from the diode on the "left" end of the log
  double m_plusEnd; // response from the diode on the "right" end of the log  

  
  unsigned int         m_index;    // index of this crystal within a layer
  unsigned int         m_layer;    // layer number
  unsigned int         m_id;	   // Log index number


};


//! Serialize the object for writing
inline StreamBuffer& IrfCalHit::serialize( StreamBuffer& s ) const                 {
  ContainedObject::serialize(s);
  return s
    << m_energy
    << m_minusEnd
    << m_plusEnd;
}


//! Serialize the object for reading
inline StreamBuffer& IrfCalHit::serialize( StreamBuffer& s )                       {
  ContainedObject::serialize(s);
  return s
    >> m_energy
    >> m_minusEnd
    >> m_plusEnd;
}


//! Fill the ASCII output stream
inline std::ostream& IrfCalHit::fillStream( std::ostream& s ) const                {
  return s
    << "class IrfCalHit :"
    << "\n    Energy    = "
    << m_energy
    << m_minusEnd
    << m_plusEnd;
}

  
//! Definition of all container types of IrfCalHit
template <class TYPE> class ObjectVector;
typedef ObjectVector<IrfCalHit>     IrfCalHitVector;
template <class TYPE> class ObjectList;
typedef ObjectList<IrfCalHit>       IrfCalHitList;


#endif // GLASTEVENT_IrfCalHit_H
