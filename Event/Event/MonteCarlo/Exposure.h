#ifndef Exposure_h
#define Exposure_h

#include <vector>
#include "GaudiKernel/ContainedObject.h"

// Include all Glast container types here
// to simplify inlude statements in algorithms
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/ObjectList.h"
#include "GaudiKernel/IInterface.h"

static const CLID& CLID_Exposure = InterfaceID("Exposure", 1, 0);

namespace Event { //Namespace
    
/** @class Exposure
* @brief TDS Container class to hold the information in the exposure (D2) Database
*
* @author Sean Robinson
*
* $Header$
*/
    class Exposure : virtual public ContainedObject
    {
    public:
        Exposure(){}
        ~Exposure(){}
        
        //! Provide ability to output some information about self...
        void writeOut(MsgStream& log) const;
        
        //! GAUDI members to be use by the converters
        static const CLID&  classID()           {return CLID_Exposure;}
        virtual const CLID& clID()        const {return classID();}
        
        //! completely initialize a new object.
        void init(
            double  intrvalstart,
            double lat, double lon,double alt,
            double posX, double  posY, double posZ,
            double RAX, double DECX,
            double RAZ, double DECZ);
        // special to initialize from objects?
        //voit init( TimeStamp time, astro::EarthCoordinate position, astro::SkiDir zenith, astro::SkyDir xaxis);
        
        //! methods to return all internal parameters.
        double RAX()const{return m_RAX;}
        double RAZ()const{return m_RAZ;}
        double DECX()const{return m_DECX;}
        double DECZ()const{return m_DECZ;}
        double lat()const{return m_lat;}
        double lon()const{return m_lon;}
        double alt()const{return m_alt;}
        double intrvalstart()const{return m_intrvalstart;}
        double posX()const{return m_posX;}
        double posY()const{return m_posY;}
        double posZ()const{return m_posZ;}

        /// Fill the ASCII output stream
        virtual std::ostream& fillStream( std::ostream& s ) const;

        
    private:
        double m_intrvalstart; /// mission elapsed time for start
        double m_lat, m_lon, m_alt; ///position  of the LAT (degrees, km)
        double m_posX,m_posY,m_posZ; /// position (in inertial celestial cartesian coordinates)
        double m_RAX, m_DECX; /// direction of LAT x-axis
        double m_RAZ, m_DECZ; /// direction of LAT z-axis
    };
    typedef ObjectVector<Exposure>  ExposureCol;
    


    inline std::ostream& Exposure::fillStream(std::ostream& log) const
    {
        log << "\tInterval Start = " << intrvalstart() 
            << "\tlat =" << lat()    << "\tlon =" << lon()   << "\talt =" << alt() 
            << "\trax =" << RAX()    << "\tdecx =" << DECX() 
            << "\traz =" << RAZ()    << "\tdecz =" << DECZ() 
            << "\tPosition = (" << posX()  << ", = " << posY()  << ", " << posZ() << ") " 
            << std::endl;
        return log;
    }

        //! completely initialize a new object.
inline void Exposure::init(
            double  intrvalstart,
            double lat, double lon,double alt,
            double posX, double  posY, double posZ,
            double RAX, double DECX,
            double RAZ, double DECZ){ 
    
    //set all the object parameters.
    m_intrvalstart = intrvalstart;
    m_lat = lat;
    m_lon = lon;
    m_alt = alt;
    m_posX = posX;
    m_posY = posY;
    m_posZ = posZ;
    m_RAX = RAX;
    m_RAZ = RAZ;
    m_DECX = DECX;
    m_DECZ = DECZ;
}
// TODO: add streaming methods
}; //Namespace

#endif