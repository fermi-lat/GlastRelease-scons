#ifndef D2Entry_h
#define D2Entry_h
/** 
* @class D2Entry
*
* @brief TDS Container class to hold the information in the D2 Database
*
* @author Sean Robinson
*
* $Header$
*/

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "Event/Recon/TkrRecon/TkrPatCand.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/ContainedObject.h"

// Include all Glast container types here
// to simplify inlude statements in algorithms
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/ObjectList.h"

/*!
//------------------------------------------------------------------------------
//
// ClassName:   D2Entry
//  
// Description: The Monte Carlo D2 database (exposure/livetime) entry structure.
//              
// Author:      Sean Robinson
//
// Changes:     [put changes here]
//------------------------------------------------------------------------------
*/


extern const CLID& CLID_D2Entry;

namespace Event { //Namespace
    
    class D2Entry : virtual public ContainedObject//public DataObject
    {
    public:
        D2Entry(){}
        ~D2Entry(){}
        
        //! Provide ability to output some information about self...
        void writeOut(MsgStream& log) const;
        
        //! GAUDI members to be use by the converters
        static const CLID&  classID()           {return CLID_D2Entry;}
        virtual const CLID& clID()        const {return classID();}
        
        //! completely initialize a new object.
        void init(double  posX,double  posY,double  posZ,double RAX,double RAZ,
            double DECX,double DECZ,double  RAZenith,double  DECZenith, double lat,
            double lon,double alt,double  intrvalstart,double  intrvalend,
            double livetime, double RAMoon, double DECMoon,double RASun,double DECSun,
            bool  SAA);
        
        //! methods to return all internal parameters.
        double RAX()const{return m_RAX;}
        double RAZ()const{return m_RAZ;}
        double DECX()const{return m_DECX;}
        double DECZ()const{return m_DECZ;}
        double lat()const{return m_lat;}
        double lon()const{return m_lon;}
        double alt()const{return m_alt;}
        double livetime()const{return m_livetime;}
        double RAMoon()const{return m_RAMoon;}
        double DECMoon()const{return m_DECMoon;}
        double RASun()const{return m_RASun;}
        double DECSun()const{return m_DECSun;}
        double intrvalstart()const{return m_intrvalstart;}
        double intrvalend()const{return m_intrvalend;}
        double posX()const{return m_posX;}
        double posY()const{return m_posY;}
        double posZ()const{return m_posZ;}
        double RAZenith()const{return m_RAZenith;}
        double DECZenith()const{return m_DECZenith;}
        bool SAA()const{return m_SAA;}
        
        ///methods to manually change a parameter.
        void RAX(double RAX){m_RAX = RAX;}
        void RAZ(double RAZ){m_RAZ = RAZ;}
        void DECX(double DECX){m_DECX = DECX;}
        void DECZ(double DECZ){m_DECZ = DECZ;}
        void lat(double lat){m_lat = lat;}
        void lon(double lon){m_lon = lon;}
        void alt(double alt){m_alt = alt;}
        void livetime(double livetime){m_livetime = livetime;}
        void RAMoon(double RAMoon){m_RAMoon = RAMoon;}
        void DECMoon(double DECMoon){m_DECMoon = DECMoon;}
        void RASun(double RASun){m_RASun = RASun;}
        void DECSun(double DECSun){m_DECSun = DECSun;}
        void intrvalstart(double intrvalstart){m_intrvalstart = intrvalstart;}
        void intrvalend(double intrvalend){m_intrvalend = intrvalend;}
        void posX(double posX){m_posX = posX;}
        void posY(double posY){m_posY = posY;}
        void posZ(double posZ){m_posZ = posZ;}
        void RAZenith(double RAZenith){m_RAZenith = RAZenith;}
        void DECZenith(double DECZenith){m_DECZenith = DECZenith;}
        void SAA(bool SAA){m_SAA = SAA;}
        
        
    private:
        double m_RAX,m_RAZ,m_DECX,m_DECZ; /// pointing characteristics of the LAT
        double m_RAZenith, m_DECZenith; /// pointing characteristic of the zenith direction.
        double m_lat,m_lon,m_alt; ///position characteristics of the LAT
        double m_livetime,m_intrvalstart,m_intrvalend; /// livetime for the current interval
        double m_RAMoon, m_DECMoon, m_RASun, m_DECSun; /// pointing of the sun and moon
        bool m_SAA; /// whether or not the satellite is in the SAA.
        double m_posX,m_posY,m_posZ; /// position (in inertial celestial cartesian coordinates)
    };
    typedef ObjectList<D2Entry>       D2EntryCol;
    
}; //Namespace

#endif