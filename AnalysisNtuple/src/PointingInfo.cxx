/** @file PointingInfo.cxx
@brief declaration and definition of the class PointingInfo


$Header$

*/
class MsgStream;
#include "Event/MonteCarlo/Exposure.h"
#include "astro/PointingInfo.h"
#include "AnalysisNtuple/PointingInfo.h"

#include "astro/SkyDir.h"
#include "astro/GPS.h"
#include "astro/PointingHistory.h"
#include "astro/JulianDate.h"
#include "astro/EarthCoordinate.h"

#include "ntupleWriterSvc/INTupleWriterSvc.h"

#include <iomanip>

using namespace AnalysisNtuple;

namespace {
    const double R2D = 180./M_PI;
	double ft2_start0 = -1.;
	double ft2_stop0  = -1.;
}

void PointingInfo::setHistoryFile( const std::string filename)
{
    m_filename = filename;
	if(filename!="") m_history = new astro::PointingHistory(m_filename);
}

void PointingInfo::execute( const astro::GPS& gps)
{
    start = gps.time();
    CLHEP::Hep3Vector pos_km = gps.position();
    CLHEP::Hep3Vector location = 1.e3* pos_km; 
    
    // cartesian location of the LAT (in m)
    sc_position[0] = location.x();
    sc_position[1] = location.y();
    sc_position[2] = location.z(); 

    ra_zenith = gps.zenithDir().ra();
    dec_zenith= gps.zenithDir().dec();

    ra_scx =    gps.xAxisDir().ra();
    dec_scx =   gps.xAxisDir().dec();

    ra_scz =    gps.zAxisDir().ra();
    dec_scz =   gps.zAxisDir().dec();

    const astro::EarthCoordinate& earth(gps.earthpos());

    lat_geo = earth.latitude(); 
    lon_geo = earth.longitude(); 
    
    rad_geo = earth.altitude(); 
    L=earth.L();
    B=earth.B();

    R=earth.R();

    lambda=earth.lambda();

    CLHEP::Hep3Vector bField = earth.magnetic_field();
    bEast  = bField.x();
    bNorth = bField.y();
    bUp    = bField.z();

    lat_mag = earth.geolat();
    in_saa= earth.insideSAA()? 1:0;
    zenith_scz = gps.zenithDir().difference(gps.zAxisDir());
 
    // sign the rocking angle: compare the declination of the zenith with that of the LAT axis
    if(dec_scz < dec_zenith) zenith_scz *= -1.0;
    zenith_scz *= R2D; // and convert to degrees

	if(m_filename!="") {
		astro::PointingHistory& history = *m_history;
		const astro::LatProperties& latProp = history(start).latProperties();
        ft2_start = latProp.interval_start();
        ft2_stop =  latProp.interval_stop();
		// no need to replace latProperties if we're in the same interval
		if(ft2_start!=ft2_start0 || ft2_stop!=ft2_stop0) {
			lat_mode   = latProp.lat_mode();
			lat_config = latProp.lat_config();
			data_qual  = latProp.data_qual();
			rock_angle = latProp.rock_angle();
	         
			livetime_frac = latProp.livetime_frac();
			ft2_start0 = ft2_start;
			ft2_stop0  = ft2_stop;

			//std::cout << "Anatup::PointingInfo - LatProperties " 
			//	<< " diff " << ft2_stop - ft2_start << " ft2 start "   
			//	<< std::setiosflags(std::ios::fixed) << std::setprecision(2) << std::setw(17) << ft2_start << " " 
			//	<< std::endl;


		}
	} else {
        lat_mode   = 0;
        lat_config = 0;
        data_qual  = 0;
        rock_angle = 0.;
         
        ft2_start = 0.;
        ft2_stop =  0.;

        livetime_frac = 0.;
	}

	//std::cout << "Anatup::PointingInfo - GPS" << std::setprecision(7) << "McIl-L " << L 
	//	<< " livetime_frac " << livetime_frac  << std::endl;
}        

/** @page anatup_vars 
    @section Pt  Pt Variables

    These items are added to the merit tuple  to give the current instrument orientation 

<table>
<tr><th> Variable <th>Type<th> Description
<tr><td> PtTime       <td>D<td> (s) Current time, same as the elapsed time
<tr><td> PtLat,PtLon  <td>F<td> (deg) lattitude and longitude
<tr><td> PtAlt        <td>F<td> (km) altitude
<tr><td> PtMagLat     <td>F<td> (deg) magnetic latitude, signed; see PtLambda
<tr><td> PtPos[3]     <td>F<td> (m) current orbit position
<tr><td> PtRax,PtDecx <td>F<td> (deg) equatorial coordinates for orientation of S/C x-axis
<tr><td> PtRaz,PtDecz <td>F<td> (deg) equatorial coordinates for orientation of S/C z-axis 
<tr><td> PtSCzenith   <td>F<td> (deg) current angle between zenith and S/C z-axis 
                                   Now signed... + means rocked north
<tr><td> PtMcIlwainB  <td>F<td> McIlwain-L parameter
<tr><td> PtMcIlwainL  <td>F<td> McIlwain-B parameter
<tr><td> PtLambda     <td>F<td> (deg)Lambda parameter, signed according to whether PlLambda 
                                increases (+) or decreases (-) with increasing PtLat
                                Test is done for a one-degree increment of PtLat
<tr><td> PtR          <td>F<td> distance to the dipole center, in Earth radii
<tr><td> PtBEast      <td>F<td> East component of the magnetic field
<tr><td> PtBNorth     <td>F<td> North component of the magnetic field
<tr><td> PtBUp        <td>F<td> Upward component of the magnetic field
<tr><td> PtLatMode    <td>I<td> LAT mode from FT2 file
<tr><td> PtLatConfig  <td>I<td> LAT configuration from FT2 file
<tr><td> PtDataQual   <td>I<td> LAT data quality word from FT2 file
<tr><td> PtRockAngle  <td>F<td> Rocking angle from FT2 file
<tr><td> PtLivetimeFrac <td>F<td> livetime fraction from FT2 file: Livetime/(Stop-Start)
</table>
*/

//------------------------------------------------------------------------
void PointingInfo::setPtTuple(INTupleWriterSvc* tuple, const std::string& tname)
{
    if( tuple==0 ) return;
    tuple->addItem(tname, "PtTime",   &start);
    tuple->addItem(tname, "PtPos[3]", sc_position);
    tuple->addItem(tname, "PtLat",    &lat_geo);
    tuple->addItem(tname, "PtLon",    &lon_geo);
    tuple->addItem(tname, "PtMagLat", &lat_mag);
    tuple->addItem(tname, "PtAlt",    &rad_geo);
    tuple->addItem(tname, "PtRaz",    &ra_scz);
    tuple->addItem(tname, "PtDecz",   &dec_scz);
    tuple->addItem(tname, "PtRax",    &ra_scx);
    tuple->addItem(tname, "PtDecx",   &dec_scx);
    tuple->addItem(tname, "PtSCzenith", &zenith_scz);
    tuple->addItem(tname, "PtMcIlwainB", &B);
    tuple->addItem(tname, "PtMcIlwainL", &L);

    tuple->addItem(tname, "PtLambda", &lambda);
    tuple->addItem(tname, "PtR",      &R);
    tuple->addItem(tname, "PtBEast",  &bEast);
    tuple->addItem(tname, "PtBNorth", &bNorth);
    tuple->addItem(tname, "PtBUp",    &bUp);

    tuple->addItem(tname, "PtLATMode",      &lat_mode);
    tuple->addItem(tname, "PtLATConfig",    &lat_config);
    tuple->addItem(tname, "PtDataQual",     &data_qual);
    tuple->addItem(tname, "PtRockAngle",    &rock_angle);
    tuple->addItem(tname, "PtLivetimeFrac", &livetime_frac);

    //tuple->addItem(tname, "PtInSAA",        &in_saa);

}
