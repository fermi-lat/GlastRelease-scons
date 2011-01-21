#ifndef Event_AcdRecon_H
#define Event_AcdRecon_H 1

#include <iostream>
#include "idents/AcdId.h"
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/StreamBuffer.h"
#include "GaudiKernel/DataObject.h"

#include "Event/TopLevel/Definitions.h"
#include "Event/Recon/AcdRecon/AcdTkrIntersection.h"
#include "Event/Recon/AcdRecon/AcdTkrPoca.h"
#include "Event/Recon/AcdRecon/AcdHit.h"
#include "Event/Recon/AcdRecon/AcdTkrHitPoca.h"
#include "Event/Recon/AcdRecon/AcdTkrGapPoca.h"
#include "Event/Recon/AcdRecon/AcdTkrPoint.h"
#include "Event/Recon/AcdRecon/AcdSplashVars.h"

#include <vector>

#include "GaudiKernel/IInterface.h"

static const CLID& CLID_AcdRecon = InterfaceID("AcdRecon", 2, 0);

/** @class Event::AcdRecon        
* @brief Reconstruction data for ACD
*
* The reconstruction data consists of:
* - Number of Acd Tiles over veto threshold.
* - Number of Acd Ribbons over veto threshold.
* - Minimum Distance of Closest Approach (DOCA)
* - List of minimum DOCAs for top and side rows.
* - Minimum Active Distance quantities (2)
* - List of minimum Active Distance quantities for top and side rows (2).
* - DOCA using the reconstructed gamma direction.
* - Collection of reconstructed energies detected by each ACD Tile.
* - Active Distance quantity for ribbons (2)
* - Total MC Energy (MeV) deposited in the whole ACD system.
* - Total MC Energy (MeV) deposited in the ribbons.
*
* - Collection of calibrated ACD signals seen in each tile
* - Collection of all track extrapolations into ACD geant model elements
* - Collection of all track extrapolations to ACD elements with signals
* - Collection of all track extrapolations to nominal location of ACD
* - Collection of all backsplash estiamtes from track extrapolations to CAL
*
*                                 
* @author Heather Kelly
* $Header$          
*/

namespace Event {
    
    class AcdRecon : virtual public DataObject  { 
        
    public:
        AcdRecon()
            : m_totEnergy(0.0),
            m_totRibbonEnergy(0.0),
            m_tileCount(0),
            m_ribbonCount(0),
            m_gammaDoca(-99999.0),
            m_doca(-99999.0)
        {};
        
        AcdRecon(double e, double ribbonEng, int count, 
            int ribbonCount, double gDoca, double doca, 
            const idents::AcdId &minDocaId, double actDist,
            const idents::AcdId &maxActDistId, 
            const std::vector<double> &rowDoca,
            const std::vector<double> &rowActDist,
            const std::vector<idents::AcdId>& idCol, 
            const std::vector<double>& energyCol,
            const std::vector<AcdTkrIntersection*>& acdTkrIntersections,
            const std::vector<AcdTkrPoca*>& acdTkrPocas,
            const std::vector<AcdHit*>& acdHits,
            const std::vector<AcdTkrHitPoca*>& acdTkrHitPocas,
            const std::vector<AcdTkrGapPoca*>& acdTkrGapPocas,
                 const std::vector<AcdTkrPoint*>& acdTkrPoints,
                 const std::vector<AcdSplashVars*>& acdSplashVars,
            double cornerDoca)
            : m_totEnergy(e),
            m_totRibbonEnergy(ribbonEng),
            m_tileCount(count),
            m_ribbonCount(ribbonCount),
            m_gammaDoca(gDoca),
            m_doca(doca),
            m_minDocaId(minDocaId),
            m_actDist(actDist),
            m_maxActDistId(maxActDistId),
            m_rowDocaCol(rowDoca),
            m_rowActDistCol(rowActDist),
            m_idCol(idCol),
            m_energyCol(energyCol),            
            m_ribbon_actDist(-2000.0),
              m_ribbon_actDist_id(idents::AcdId(0,0)),            
            m_acdTkrIntersections(acdTkrIntersections),
            m_acdTkrPocas(acdTkrPocas),
               m_acdHits(acdHits),
                m_cornerDoca(cornerDoca),
            m_acdTkrHitPocas(acdTkrHitPocas),
            m_acdTkrGapPocas(acdTkrGapPocas),          
            m_acdTkrPoints(acdTkrPoints),
            m_acdSplashVars(acdSplashVars)
          {
            m_actDist3D = -2000.0;
            m_rowActDist3DCol.resize(4, -2000.0);
        };

        // DC: this was, as compared to the previous one,
        // has two additionnal ribbon* arguments, and
        // three 
        AcdRecon(double e, double ribbonE, int count, int ribbonCount, 
            double gDoca, double doca, 
            const idents::AcdId &minDocaId, double actDist,
            const idents::AcdId &maxActDistId, 
            const std::vector<double> &rowDoca,
            const std::vector<double> &rowActDist,
            const std::vector<idents::AcdId>& idCol, 
            const std::vector<double>& energyCol,
            double ribbon_actDist, const idents::AcdId ribbon_actDist_id,
            const std::vector<AcdTkrIntersection*>& acdTkrIntersections,
            const std::vector<AcdTkrPoca*>& acdTkrPocas,
                 const std::vector<AcdHit*>& acdHits,
                 const std::vector<AcdTkrHitPoca*>& acdTkrHitPocas,
            const std::vector<AcdTkrGapPoca*>& acdTkrGapPocas,
                 const std::vector<AcdTkrPoint*>& acdTkrPoints,
                 const std::vector<AcdSplashVars*>& acdSplashVars,
            double actDist3D, const idents::AcdId &maxActDist3DId, 
            const std::vector<double> &rowActDist3D, double cornerDoca)
            : m_totEnergy(e),
            m_totRibbonEnergy(ribbonE),
            m_tileCount(count),
            m_ribbonCount(ribbonCount),
            m_gammaDoca(gDoca),
            m_doca(doca),
            m_minDocaId(minDocaId),
            m_actDist(actDist),
            m_actDist3D(actDist3D),
            m_actDist3D_down(0.0),
            m_maxActDistId(maxActDistId),
            m_maxActDist3DId(maxActDist3DId),
            m_rowDocaCol(rowDoca),
            m_rowActDistCol(rowActDist),
            m_rowActDist3DCol(rowActDist3D),
            m_idCol(idCol),
            m_energyCol(energyCol),            
            m_ribbon_actDist(ribbon_actDist),
            m_ribbon_actDist_id(ribbon_actDist_id),
            m_acdTkrIntersections(acdTkrIntersections),
            m_acdTkrPocas(acdTkrPocas),
            m_acdHits(acdHits),
            m_cornerDoca(cornerDoca),
            m_acdTkrHitPocas(acdTkrHitPocas),
            m_acdTkrGapPocas(acdTkrGapPocas),
            m_acdTkrPoints(acdTkrPoints),
            m_acdSplashVars(acdSplashVars)          
            {
                m_rowActDist3DCol_down.clear();

            };


        
        virtual ~AcdRecon() { };


        void init(double e, double ribbonE, int count, 
            int ribbonCount, double gDoca, double doca, 
            const idents::AcdId &minDocaId, double actDist,
            const idents::AcdId &maxActDistId, 
            const std::vector<double> &rowDoca,
            const std::vector<double>&  rowActDist,
            const std::vector<idents::AcdId>& idCol, 
            const std::vector<double>&  energyCol,
            double actDist3D, const idents::AcdId &maxActDist3DId,
            const std::vector<double> &rowActDist3D,
            const std::vector<AcdTkrIntersection*>& acdTkrIntersections,
                  const std::vector<AcdTkrPoca*>& acdTkrPocas,
                  const std::vector<AcdHit*>& acdHits,
                  const std::vector<AcdTkrHitPoca*>& acdTkrHitPocas,
                  const std::vector<AcdTkrGapPoca*>& acdTkrGapPocas,
                  const std::vector<AcdTkrPoint*>& acdTkrPoints,
                  const std::vector<AcdSplashVars*>& acdSplashVars,
            double ribbonActDist=2000.0, 
            const idents::AcdId &ribActDistId=idents::AcdId(0,0),
            double cornerDoca=2000.0) ;

        void initActDist3D_down(double actDist3D,
                                const idents::AcdId &maxActDist3DId,
                                const std::vector<double> &rowActDist3D) {

            m_actDist3D_down = actDist3D;
            m_maxActDist3DId_down = maxActDist3DId;
            m_rowActDist3DCol_down = rowActDist3D;
        }

        void clear();

        //! Retrieve reference to class definition structure
        virtual const CLID& clID() const    { return AcdRecon::classID(); }
        static const  CLID& classID()       { return CLID_AcdRecon; }
        
        inline const double getEnergy()     const { return m_totEnergy; };
        inline const double getRibbonEnergy() const { return m_totRibbonEnergy; };
        inline const int    getTileCount()  const { return m_tileCount; };
        inline const int    getRibbonCount()  const { return m_ribbonCount; };
        inline const double getGammaDoca()  const { return m_gammaDoca; };
        inline const double getDoca()       const { return m_doca; };
        inline const double getCornerDoca() const { return m_cornerDoca; };
        inline const double getActiveDist() const { return m_actDist; };
        inline const double getActiveDist3D() const { return m_actDist3D; };
        inline const double getActiveDist3D_Up() const { return m_actDist3D; };
        inline const double getActiveDist3D_Down() const { return m_actDist3D_down; };
        inline const double getRibbonActiveDist()            const { return m_ribbon_actDist; };
        inline const idents::AcdId& getRibbonActiveDistId()  const { return m_ribbon_actDist_id; };
        inline const idents::AcdId& getMinDocaId()           const { return m_minDocaId; };
        inline const idents::AcdId& getMaxActDistId()        const { return m_maxActDistId; };
        inline const idents::AcdId& getMaxActDist3DId()        const { return m_maxActDist3DId; };
        inline const idents::AcdId& getMaxActDist3DId_Up()        const { return m_maxActDist3DId; };
        inline const idents::AcdId& getMaxActDist3DId_Down()        const { return m_maxActDist3DId_down; };
        inline const std::vector<double>& getRowDocaCol()    const { return m_rowDocaCol; };
        inline const std::vector<double>& getRowActDistCol() const { return m_rowActDistCol; };
        inline const std::vector<double>& getRowActDist3DCol() const { return m_rowActDist3DCol; };
        inline const std::vector<double>& getRowActDist3DCol_Up() const { return m_rowActDist3DCol; };
        inline const std::vector<double>& getRowActDist3DCol_Down() const { return m_rowActDist3DCol_down; };
        inline const std::vector<idents::AcdId>& getIdCol()  const { return m_idCol; };
        inline const std::vector<double>& getEnergyCol()     const { return m_energyCol; };
        inline const AcdTkrIntersectionCol& getAcdTkrIntersectionCol() const { return m_acdTkrIntersections; };
        inline const AcdTkrPocaCol& getAcdTkrPocaCol() const { return m_acdTkrPocas; };
        inline const AcdHitCol& getAcdHitCol() const { return m_acdHits; };
        inline const AcdTkrHitPocaCol& getAcdTkrHitPocaCol() const { return m_acdTkrHitPocas; };
        inline const AcdTkrGapPocaCol& getAcdTkrGapPocaCol() const { return m_acdTkrGapPocas; };
        inline const AcdTkrPointCol& getAcdTkrPointCol() const { return m_acdTkrPoints; }
        inline const AcdSplashVarsCol& getAcdSplashVarsCol() const { return m_acdSplashVars; }

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
        
        /// Total MC energy in MeV deposited in the whole ACD system
        /// This remains as a check on MC runs
        double m_totEnergy;
        /// Total MeV energy deposited in the ribbons
        double m_totRibbonEnergy;
        /// Total number of ACD tiles above threshold
        int    m_tileCount;
        /// Total number of ACD ribbons above threshold 
        int m_ribbonCount;

        /// Distance of Closest Approach for the reconstructed gamma, if there is one
        double m_gammaDoca;

        double m_doca;

        /// record of the tile with the minimum Distance of Closest Approach
        idents::AcdId m_minDocaId;

        /// DOCA calculation using edge of tiles (Bill Atwood)
        double m_actDist;
        /// Active Distance calc for both upward and downward tracks
        double m_actDist3D, m_actDist3D_down;
    
        /// record of the tile with the maximum Active Distance 
        idents::AcdId m_maxActDistId;
        idents::AcdId m_maxActDist3DId, m_maxActDist3DId_down;
 

        /// Collection of distance of closest approach calculations
        /// for each side row of the ACD
        ///   zeroth element corresponds to the top, and 
        ///   index 1 corresponds to the first row closest to the top, etc.
        std::vector<double>  m_rowDocaCol;
        /// Collection of Active Distance quantities for each row of the ACD
        ///    zeroth element corresponds to the top, and 
        ///    index 1 corresponds to the first row closest to the top, etc.
        std::vector<double>  m_rowActDistCol;
        std::vector<double>  m_rowActDist3DCol;
        std::vector<double>  m_rowActDist3DCol_down;
        
        /// Reconstructed energy per ACD digi - MC for now
        std::vector<idents::AcdId> m_idCol;
        std::vector<double>        m_energyCol;

        /// Active Distance calculation for ribbons
        double         m_ribbon_actDist;
        /// Id of the ribbon corresponding to the Active Distance
        idents::AcdId  m_ribbon_actDist_id;        

        /// the vector of track intersections w/ the acd
        AcdTkrIntersectionCol m_acdTkrIntersections;

        /// the vector of track poca w/ the acd
        AcdTkrPocaCol m_acdTkrPocas;        

        /// the vector of calibrated ACD hits
        AcdHitCol m_acdHits;        

        /// Bill's variable to measure DOCA to rays along corner side gaps
        double m_cornerDoca;

        /// new calculation of tkr pocas
        AcdTkrHitPocaCol m_acdTkrHitPocas;

        /// pocas to the gaps
        AcdTkrGapPocaCol m_acdTkrGapPocas;

        /// points where the track cross the nomial ACD
        AcdTkrPointCol m_acdTkrPoints;

        /// Store the angles and such to study backsplash
        AcdSplashVarsCol m_acdSplashVars;


    };


    inline void AcdRecon::clear() {
        m_rowDocaCol.clear();
        m_rowActDistCol.clear();
        m_rowActDist3DCol.clear();
        m_rowActDist3DCol_down.clear();
        m_idCol.clear();
        m_energyCol.clear();
        m_acdTkrIntersections.clear();
        m_acdTkrPocas.clear();
        m_acdHits.clear();
        m_acdTkrHitPocas.clear();
        m_acdTkrGapPocas.clear();
        m_acdTkrPoints.clear();
        m_acdSplashVars.clear();
    }


    inline void AcdRecon::init(double e, double ribbonE, int count, 
            int ribbonCount, double gDoca, double doca, 
            const idents::AcdId &minDocaId, double actDist,
            const idents::AcdId &maxActDistId, 
            const std::vector<double> &rowDoca,
            const std::vector<double>&  rowActDist,
            const std::vector<idents::AcdId>& idCol, 
            const std::vector<double>& energyCol, 
            double actDist3D, const idents::AcdId &maxActDist3DId,
            const std::vector<double> &rowActDist3D,
            const std::vector<AcdTkrIntersection*>& acdTkrIntersections,
            const std::vector<AcdTkrPoca*>& acdTkrPocas,
                               const std::vector<AcdHit*>& acdHits,
                               const std::vector<AcdTkrHitPoca*>& acdTkrHitPocas,
                               const std::vector<AcdTkrGapPoca*>& acdTkrGapPocas,    
                               const std::vector<AcdTkrPoint*>& acdTkrPoints,
                               const std::vector<AcdSplashVars*>& acdSplashVars,
            double ribbon_actDist, const idents::AcdId &ribbonId, 
            double cornerDoca)
    {
        m_totEnergy  = e;
        m_totRibbonEnergy = ribbonE;
        m_tileCount  = count;
        m_ribbonCount  = ribbonCount;
        m_gammaDoca  = gDoca;
        m_doca       = doca;
        m_actDist    = actDist;
        m_actDist3D    = actDist3D;
        m_actDist3D_down = 0.0;
        m_ribbon_actDist    = ribbon_actDist;
        m_ribbon_actDist_id = ribbonId;
        m_minDocaId  = minDocaId;
        m_maxActDistId  = maxActDistId;
        m_maxActDist3DId  = maxActDist3DId;
        m_rowDocaCol = rowDoca;
        m_rowActDistCol = rowActDist;
        m_rowActDist3DCol = rowActDist3D;
        m_rowActDist3DCol_down.clear();
        m_idCol      = idCol;
        m_energyCol  = energyCol;
        m_acdTkrIntersections.clear();
        for ( std::vector<AcdTkrIntersection*>::const_iterator itr = acdTkrIntersections.begin();
              itr != acdTkrIntersections.end(); itr++ ) {
          AcdTkrIntersection* iSect = const_cast<AcdTkrIntersection*>(*itr);
          m_acdTkrIntersections.add(iSect);
        }
        m_acdTkrPocas.clear();
        for ( std::vector<AcdTkrPoca*>::const_iterator itrPoca = acdTkrPocas.begin();
              itrPoca != acdTkrPocas.end(); itrPoca++ ) {
          AcdTkrPoca* poca = const_cast<AcdTkrPoca*>(*itrPoca);
          m_acdTkrPocas.add(poca);
        }
        m_acdHits.clear();
        for ( std::vector<AcdHit*>::const_iterator itrHit = acdHits.begin();
              itrHit != acdHits.end(); itrHit++ ) {
          AcdHit* hit = const_cast<AcdHit*>(*itrHit);
          m_acdHits.add(hit);
        }        
        m_acdTkrHitPocas.clear();
        for ( std::vector<AcdTkrHitPoca*>::const_iterator itrHitPoca = acdTkrHitPocas.begin();
              itrHitPoca != acdTkrHitPocas.end(); itrHitPoca++ ) {
          AcdTkrHitPoca* hitPoca = const_cast<AcdTkrHitPoca*>(*itrHitPoca);
          m_acdTkrHitPocas.add(hitPoca);
        }
        m_acdTkrGapPocas.clear();
        for ( std::vector<AcdTkrGapPoca*>::const_iterator itrGapPoca = acdTkrGapPocas.begin();
              itrGapPoca != acdTkrGapPocas.end(); itrGapPoca++ ) {
          AcdTkrGapPoca* gapPoca = const_cast<AcdTkrGapPoca*>(*itrGapPoca);
          m_acdTkrGapPocas.add(gapPoca);
        }
        m_acdTkrPoints.clear();
        for ( std::vector<AcdTkrPoint*>::const_iterator itrPoint = acdTkrPoints.begin();
              itrPoint != acdTkrPoints.end(); itrPoint++ ) {
          AcdTkrPoint* point = const_cast<AcdTkrPoint*>(*itrPoint);
          m_acdTkrPoints.add(point);
        }
        m_acdSplashVars.clear();
        for ( std::vector<AcdSplashVars*>::const_iterator itrSplash = acdSplashVars.begin();
              itrSplash != acdSplashVars.end(); itrSplash++ ) {
          AcdSplashVars* splash = const_cast<AcdSplashVars*>(*itrSplash);
          m_acdSplashVars.add(splash);
        }
        m_cornerDoca = cornerDoca;
    }

    
    /// Serialize the object for writing
    inline StreamBuffer& AcdRecon::serialize( StreamBuffer& s ) const
    {
        DataObject::serialize(s);
        return s
            << m_totEnergy
            << m_totRibbonEnergy
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
