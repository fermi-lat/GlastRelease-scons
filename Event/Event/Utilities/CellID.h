// $Header$
#ifndef LHCBEVENT_CELLID_H
#define LHCBEVENT_CELLID_H 1


// Include files
#include <iostream>
#include "GaudiKernel/StreamBuffer.h"
#include "GlastEvent/TopLevel/Definitions.h"


//------------------------------------------------------------------------------
//
// ClassName:   CellID
//  
// Description: Cell identifier of the calorimeter cells
//
// Author:      Pavel Binko
// Changes:     P.Binko 19/10/1999 : Formating of ASCII output
//
//------------------------------------------------------------------------------


class CellID                                                                   {

  // Encoding taken from SICB
  // Restriction : CalorimeterEncoding > RawEncoding > ColumnEncoding
  //               If not, the code has to be changed
  enum { CalorimeterEncoding=1000000, RawEncoding=1000, ColumnEncoding=1 };

public:

  /// Constructors
  CellID()
    : m_id(0)                                                                { }
  CellID( long cellID )
    : m_id(cellID)                                                           { }
  CellID( long calorimeterNumber, long rawNumber, long columnNumber )          {
    setID( calorimeterNumber, rawNumber, columnNumber );
  }
  /// Destructor
  ~CellID()                                                                  { }

  /// Retrieve calorimeter number
  long calorimeterNumber() const                                               {
    return long( m_id / CalorimeterEncoding );
  }
  /// Retrieve raw number
  long rawNumber() const                                                       {
    return long( (m_id - calorimeterNumber()*CalorimeterEncoding)
      / RawEncoding );
  }
  /// Retrieve column number
  long columnNumber() const                                                    {
    return long( (m_id - calorimeterNumber()*CalorimeterEncoding
      - rawNumber()*RawEncoding)
      / ColumnEncoding );
  }
  /// Retrieve calorimeter cell ID
  long id() const                                                              {
    return m_id;
  }
  /// Update calorimeter cell ID
  void setID( long calorimeterNumber, long rawNumber, long columnNumber )      {
    m_id  = CalorimeterEncoding*calorimeterNumber
      + RawEncoding*rawNumber
      + ColumnEncoding*columnNumber;
  }

  /// Serialize the object for writing
  friend StreamBuffer& operator<< ( StreamBuffer& s, const CellID& obj )       {
    return s << obj.m_id;
  }
  /// Serialize the object for reading
  friend StreamBuffer& operator>> ( StreamBuffer& s, CellID& obj )             {
    return s >> obj.m_id;
  }

  /// Output operator (ASCII)
  friend std::ostream& operator<< ( std::ostream& s, const CellID& obj )       {
    return obj.fillStream(s);
  }
  /// Fill the output stream (ASCII)
    std::ostream& fillStream( std::ostream& s ) const                          {
    return s
      << "class CellID (calorimeter, raw, column) : ( "
      << GlastEventField( GlastEvent::field5 )
      << calorimeterNumber()
      << ", "
      << GlastEventField( GlastEvent::field5 )
      << rawNumber()
      << ", "
      << GlastEventField( GlastEvent::field5 )
      << columnNumber() << " )";
  }

private:

  // Cell identifier
  long m_id;

};


#endif    // LHCBEVENT_CELLID_H
