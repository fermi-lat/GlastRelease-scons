#ifndef GlastEvent_AcdRecon_H
#define GlastEvent_AcdRecon_H 1


// Include files
#include <iostream>
#include "idents/AcdId.h"
#include "data/IVetoData.h"
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/StreamBuffer.h"
#include "GaudiKernel/DataObject.h"

#include "GlastEvent/TopLevel/Definitions.h"

#include <vector>
#include <map>

/*!
//------------------------------------------------------------------------------
//
\class   AcdRecon        
  
\brief Reconstruction data for ACD                                
              

//------------------------------------------------------------------------------
 */

extern const CLID& CLID_AcdRecon;


class AcdRecon : virtual public DataObject  { 

public:
    AcdRecon()
        : m_totEnergy(0.0),
        m_tileCount(0),
        m_gammaDOCA(-99999.0),
        m_DOCA(-99999.0)
    {};

    AcdRecon(double e, int count, double gDoca, double doca,
        IVetoData::Tile hitTile, std::vector<double> &rowDOCA, 
        std::map<idents::AcdId,double> &energies)       
        : m_totEnergy(e),
        m_tileCount(count),
        m_gammaDOCA(gDoca),
        m_DOCA(doca),
        m_hitTile(hitTile),
        m_rowDOCA_vec(rowDOCA),
        m_energies(energies)

    {};

    /// Destructor
    virtual ~AcdRecon() { };

    //! Retrieve reference to class definition structure
    virtual const CLID& clID() const   { return AcdRecon::classID(); }
    static const CLID& classID()       { return CLID_AcdRecon; }

    inline const double energy() const { return m_totEnergy; };
    inline const int tileCount() const { return m_tileCount; };
    inline const double gammaDOCA() const { return m_gammaDOCA; };
    inline const double DOCA() const { return m_DOCA; };



    /// Serialize the object for writing
    virtual StreamBuffer& serialize( StreamBuffer& s ) const;
    /// Serialize the object for reading
    virtual StreamBuffer& serialize( StreamBuffer& s );
    /// Fill the ASCII output stream
    virtual std::ostream& fillStream( std::ostream& s ) const;



private:

    double m_totEnergy;
    int m_tileCount;
    double m_gammaDOCA;
    double m_DOCA;
    std::vector<double> m_rowDOCA_vec;

    // record of the tile with the minimum Distance of Closest Approach
    IVetoData::Tile m_hitTile;

    std::map<idents::AcdId,double> m_energies;

};


/// Serialize the object for writing
inline StreamBuffer& AcdRecon::serialize( StreamBuffer& s ) const
{
  DataObject::serialize(s);
  return s
    << m_totEnergy
    << m_tileCount
    << m_gammaDOCA
    << m_DOCA;
}


/// Serialize the object for reading
inline StreamBuffer& AcdRecon::serialize( StreamBuffer& s )
{
  DataObject::serialize(s);
  
  s >> m_totEnergy
    >> m_tileCount
    >> m_gammaDOCA
    >> m_DOCA;

  return s;
}


/// Fill the ASCII output stream
inline std::ostream& AcdRecon::fillStream( std::ostream& s ) const
{
  return s
    << "    base class AcdRecon :"
    << "\n        total energy      = "
    << EventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_totEnergy << ", "
    << "\n        tile Count              = "
    << EventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_tileCount   << " )"
    << "\n        gamma DOCA     = "
    << m_gammaDOCA << " )"
    << "\n        DOCA     = "
    << EventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_DOCA << " )";
}




#endif    // GlastEvent_AcdRecon_H

