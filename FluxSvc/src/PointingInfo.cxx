/** @file PointingInfo.cxx
@brief declaration and definition of the class PointingInfo


$Header$

*/
#include "FluxSvc/PointingInfo.h"
#include "astro/GPS.h"
#include "astro/SkyDir.h"
#include "astro/JulianDate.h"
#include "astro/EarthCoordinate.h"
#include "ntupleWriterSvc/INTupleWriterSvc.h"

#include "Event/MonteCarlo/Exposure.h"


void PointingInfo::set(double time)
{
    // Purpose: save the current status until the next tick
    using namespace astro;
 
    start = time;
    // The GPS singleton has current time and orientation
    GPS* gps = GPS::instance();
    gps->getPointingCharacteristics(time); // sets time for other functions
    Hep3Vector pos_km = gps->position(time);
    Hep3Vector location = 1.e3* pos_km; // special, needs its own time
    
    // cartesian location of the LAT (in m)
    sc_position[0] = location.x();
    sc_position[1] = location.y();
    sc_position[2] = location.z(); 

    ra_zenith = gps->RAZenith();
    dec_zenith= gps->DECZenith();
    ra_scx =    gps->RAX();
    dec_scx =   gps->DECX();

    ra_scz =    gps->RAZ();
    dec_scz =   gps->DECZ();
    //uncomment for debug check double check=astro::SkyDir(rax, decx)().dot(astro::SkyDir(raz, decz)());

    lat_geo = gps->lat(); 
    lon_geo = gps->lon(); 

    // override altitude by using shape of earth.
    EarthCoordinate loc(pos_km, JulianDate(time));
    rad_geo = loc.altitude(); 


}
Event::Exposure* PointingInfo::forTDS()const
{

    Event::Exposure* entry = new Event::Exposure;
    entry->init(start,lat_geo,lon_geo,rad_geo,
        sc_position[0],sc_position[1],sc_position[2],ra_scx,dec_scx,ra_scz,dec_scz);
    return entry;

}
//------------------------------------------------------------------------
void PointingInfo::finish(double stop_time, double live)
{
    stop = stop_time;
    livetime = live;
}
//------------------------------------------------------------------------
void PointingInfo::setFT2Tuple(INTupleWriterSvc* tuple, const std::string& tname)
{
    if( tuple==0 ) return;
    tuple->addItem(tname, "start",     &start);
    tuple->addItem(tname, "stop",      &stop);
    tuple->addItem(tname, "sc_position[3]", sc_position);
    tuple->addItem(tname, "lat_geo",   &lat_geo);
    tuple->addItem(tname, "lon_geo",   &lon_geo);
    tuple->addItem(tname, "rad_geo",   &rad_geo);
    tuple->addItem(tname, "ra_zenith", &ra_zenith);
    tuple->addItem(tname, "dec_zenith",&dec_zenith);
    tuple->addItem(tname, "ra_scz",    &ra_scz);
    tuple->addItem(tname, "dec_scz",   &dec_scz);
    tuple->addItem(tname, "ra_scx",    &ra_scx);
    tuple->addItem(tname, "dec_scx",   &dec_scx);
    tuple->addItem(tname, "livetime",  &livetime);
}
//------------------------------------------------------------------------
void PointingInfo::setPtTuple(INTupleWriterSvc* tuple, const std::string& tname)
{

       /** @page point_info pointing information
     These items are added to the merit tuple in MeritAlg, to give the current instrument orientation 


    - PtTime (s) Current time, same as the elapsed time
    - PtLat,PtLon (deg) lattitude and longitude
    - PtAlt (km) altitude
    - PtPos[3] (m) current orbit position
    - PtRax,PtDecx (deg) equatorial coordinates for orientation of S/C x-axis
    - PtRaz,PtDecz (deg) equatorial coordinates for orientation of S/C z-axis 
    */

    if( tuple==0 ) return;
    tuple->addItem(tname, "PtTime",   &start);
    tuple->addItem(tname, "PtPos[3]", sc_position);
    tuple->addItem(tname, "PtLat",    &lat_geo);
    tuple->addItem(tname, "PtLon",    &lon_geo);
    tuple->addItem(tname, "PtAlt",    &rad_geo);
    tuple->addItem(tname, "PtRaz",    &ra_scz);
    tuple->addItem(tname, "PtDecz",   &dec_scz);
    tuple->addItem(tname, "PtRax",    &ra_scx);
    tuple->addItem(tname, "PtDecx",   &dec_scx);

}