#ifndef Event_AcdRecon_H
#define Event_AcdRecon_H 1

#include <iostream>
#include "idents/AcdId.h"
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/StreamBuffer.h"
#include "GaudiKernel/DataObject.h"

#include "Event/TopLevel/Definitions.h"

#include <vector>

extern const CLID& CLID_AcdRecon;

/** @class AcdRecon        
* @brief Reconstruction data for ACD
*
* The reconstruction data consists of:
* - Total Energy deposited in the whole ACD system.
* - Number of Acd Tiles over veto threshold.
* - Minimum Distance of Closest Approach (DOCA)
* - List of minimum DOCAs for top and side rows.
* - Minimum Active Distance quantity
* - List of minimum Active Distance quantities for top and side rows.
* - DOCA using the reconstructed gamma direction.
* - Collection of reconstructed energies detected by each ACD Tile.
* - Active Distance quantity for ribbons
*                                 
* @author Heather Kelly
* $Header$          
*/
namespace Event {
    
    class AcdRecon : virtual public DataObject  { 
        
    public:
        AcdRecon()
            : m_totEnergy(0.0),
            m_tileCount(0),
            m_gammaDoca(-99999.0),
            m_doca(-99999.0)
        {};
        
        AcdRecon(double e, int count, double gDoca, double doca, double actDist,
            const idents::AcdId &minDocaId, const std::vector<double> &rowDoca,
            const std::vector<double> &rowActDist,
            const std::vector<idents::AcdId>& idCol, 
            const std::vector<double>& energyCol)
            : m_totEnergy(e),
            m_tileCount(count),
            m_gammaDoca(gDoca),
            m_doca(doca),
            m_actDist(actDist),
            m_minDocaId(minDocaId),
            m_rowDocaCol(rowDoca),
            m_rowActDistCol(rowActDist),
			m_idCol(idCol),
            m_energyCol(energyCol),
            m_ribbon_actDist(-2000.0),
            m_ribbon_actDist_id(idents::AcdId(0,0))
            
        {};


        AcdRecon(double e, int count, double gDoca, double doca, double actDist,
            const idents::AcdId &minDocaId, const std::vector<double> &rowDoca,
            const std::vector<double> &rowActDist,
            const std::vector<idents::AcdId>& idCol, 
            const std::vector<double>& energyCol,
            double ribbon_actDist, const idents::AcdId ribbon_actDist_id)
            : m_totEnergy(e),
            m_tileCount(count),
            m_gammaDoca(gDoca),
            m_doca(doca),
            m_actDist(actDist),
            m_minDocaId(minDocaId),
            m_rowDocaCol(rowDoca),
            m_rowActDistCol(rowActDist),
			m_idCol(idCol),
            m_energyCol(energyCol),
            m_ribbon_actDist(ribbon_actDist),
            m_ribbon_actDist_id(ribbon_actDist_id)
            
        {};


        
        virtual ~AcdRecon() { };

        void init(double e, int count, double gDoca, double doca, double actDist,
            const idents::AcdId &minDocaId, const std::vector<double> &rowDoca,
            const std::vector<double> &rowActDist,
            const std::vector<idents::AcdId>& idCol, 
            const std::vector<double>& energyCol,
            double ribbonActDist=2000.0, const idents::AcdId &ribActDistId=idents::AcdId(0,0));

        void clear();

        //! Retrieve reference to class definition structure
        virtual const CLID& clID() const   { return AcdRecon::classID(); }
        static const CLID& classID()       { return CLID_AcdRecon; }
        
        inline const double getEnergy() const { return m_totEnergy; };
        inline const int getTileCount() const { return m_tileCount; };
        inline const double getGammaDoca() const { return m_gammaDoca; };
        inline const double getDoca() const { return m_doca; };
        inline const double getActiveDist() const { return m_actDist; };
        inline const double getRibbonActiveDist() const { return m_ribbon_actDist; };
        inline const idents::AcdId& getRibbonActiveDistId() const { return m_ribbon_actDist_id; };
        inline const idents::AcdId& getMinDocaId() const { return m_minDocaId; };
        inline const std::vector<double>& getRowDocaCol() const { return m_rowDocaCol; };
        inline const std::vector<double>& getRowActDistCol() const { return m_rowActDistCol; };
		inline const std::vector<idents::AcdId>& getIdCol() const { return m_idCol; };
        inline const std::vector<double>& getEnergyCol() const { return m_energyCol; };

        /// Serialize the object for writing
        virtual StreamBuffer& serialize( StreamBuffer& s ) const;
        /// Serialize the object for reading
        virtual StreamBuffer& serialize( StreamBuffer& s );
        
        
        friend std::ostream& operator << (std::ostream& s, const AcdRecon& obj)
        {
            return obj.fillStream(s);
        };
        
        /// Fill the ASCII output stream
        virtual std::ostream& fillStream( std::ostream& s ) const;
        
    private:
        
        /// Total energy in MeV deposited in the whole ACD system
        double m_totEnergy;
        /// Total number of ACD tiles above threshold
        int m_tileCount;
        /// Distance of Closest Approach for the reconstructed gamma, 
        /// if there is one
        double m_gammaDoca;
        /// Minimum Distance of Closest Approach for all tracks and all ACD tiles
        double m_doca;
        /// New Bill Atwood DOCA calculation using edge of tiles
        double m_actDist;
        // Active Distance calculation for ribbons
        double m_ribbon_actDist;
        // Id of the ribbon corresponding to the Active Distance
        idents::AcdId m_ribbon_actDist_id;
        // record of the tile with the minimum Distance of Closest Approach
        idents::AcdId m_minDocaId;
        /// Collection of distance of closest approach calculations
        /// for each side row of the ACD
        /// zeroth element corresponds to the top, and index one corresponds to
        /// the first row closest to the top, etc.
        std::vector<double> m_rowDocaCol;
        /// Collection of Active Distance quantities for each row of the ACD
        /// zeroth element corresponds to the top, and index one corresponds to
        /// the first row closest to the top, etc.
        std::vector<double> m_rowActDistCol;
        
        /// Stores reconstructed energy per ACD digi
        //std::map<idents::AcdId, double> m_energyCol;
        std::vector<idents::AcdId> m_idCol;
        std::vector<double> m_energyCol;
        
    };

    inline void AcdRecon::clear() {
        m_rowDocaCol.clear();
        m_rowActDistCol.clear();
        m_idCol.clear();
        m_energyCol.clear();
    }

    inline void AcdRecon::init(double e, int count, double gDoca, double doca, double actDist,
            const idents::AcdId &minDocaId, const std::vector<double> &rowDoca,
            const std::vector<double> &rowActDist,
            const std::vector<idents::AcdId>& idCol, 
            const std::vector<double>& energyCol, double ribbon_actDist, const idents::AcdId &ribbonId)
    {
        m_totEnergy = e;
        m_tileCount = count;
        m_gammaDoca = gDoca;
        m_doca = doca;
        m_actDist = actDist;
        m_ribbon_actDist = ribbon_actDist;
        m_ribbon_actDist_id = ribbonId;
        m_minDocaId = minDocaId;
        m_rowDocaCol = rowDoca;
        m_rowActDistCol = rowActDist;
        m_idCol = idCol;
        m_energyCol = energyCol;
    }

    
    /// Serialize the object for writing
    inline StreamBuffer& AcdRecon::serialize( StreamBuffer& s ) const
    {
        DataObject::serialize(s);
        return s
            << m_totEnergy
            << m_tileCount
            << m_gammaDoca
            << m_doca;
    }
    
    
    /// Serialize the object for reading
    inline StreamBuffer& AcdRecon::serialize( StreamBuffer& s )
    {
        DataObject::serialize(s);
        
        s >> m_totEnergy
            >> m_tileCount
            >> m_gammaDoca
            >> m_doca;
        
        return s;
    }
    
    
    /// Fill the ASCII output stream
    inline std::ostream& AcdRecon::fillStream( std::ostream& s ) const
    {
        return s
            << "    base class AcdRecon :"
            << "\n        total energy      = "
            << EventFloatFormat( EventFormat::width, EventFormat::precision )
            << m_totEnergy << ", "
            << "\n        tile Count              = "
            << EventFloatFormat( EventFormat::width, EventFormat::precision )
            << m_tileCount   << " )"
            << "\n        gamma DOCA     = "
            << m_gammaDoca << " )"
            << "\n        DOCA     = "
            << EventFloatFormat( EventFormat::width, EventFormat::precision )
            << m_doca << " )"
            << "\n        ribbon Active Distance = "
            << EventFloatFormat( EventFormat::width, EventFormat::precision )
            << m_ribbon_actDist << " )";
    }
    
    
} // namespace Event

#endif    // Event_AcdRecon_H

