// $Header$

// Include files
#include <iostream>
//#include "geometry/Hep3Vector.h"
#include "Event/Utilities/SkyDir.h"
#include "geometry/CoordTransform.h"

SkyDir::SkyDir(double param1, double param2, coordSystem inputType){
    if(inputType == GALACTIC){
       double  m_l = param1;
       double  m_b = param2;
        
        if (m_l==0.){m_l+=0.000000000001;}  //to fix divide-by-zero errors
        
        //here we construct the cartesian galactic vector
        Hep3Vector gamgal(sin(m_l*M_2PI/360.)*cos(m_b*M_2PI/360.) , sin(m_b * M_2PI/360.) , cos(m_l*M_2PI/360.)*cos(m_b*M_2PI/360.));
        
        //get the transformation matrix from galactic to celestial coordinates.
        Rotation galToCel=celToGal().inverse();
        
        //and do the transform to get the cartesian celestial vector
        m_dir = galToCel*gamgal;
        
    }else if(inputType == CELESTIAL){
        double m_ra = param1;
        double m_dec = param2;
        
        if (m_ra==0.){m_ra+=0.000000000001;}  //to fix divide-by-zero errors
        
        //here we construct the cartesian celestial vector
        m_dir = Hep3Vector(sin(m_ra*M_2PI/360.)*cos(m_dec*M_2PI/360.) , sin(m_dec*M_2PI/360.) , cos(m_ra*M_2PI/360.)*cos(m_dec*M_2PI/360.));        
        
    }else{
        //improper coordinate system declaration - default things and say so.
        std::cout << "Improper coordinate System declaration in SkyDir" << std::endl;
        
        m_dir = Hep3Vector(0,0,1);
    }
}

SkyDir::SkyDir(Hep3Vector dir):
m_dir(dir){
}

Rotation SkyDir::celToGal()const{
    //gal is the matrix which rotates (cartesian)celestial coordiantes into (cartesian)galactic ones
    Rotation gal;
    double degsPerRad = 180./M_PI;
    gal.rotateZ(-282.25/degsPerRad).rotateX(-62.6/degsPerRad).rotateZ(33./degsPerRad);
    return gal; 
}


std::pair<double,double> SkyDir::setGalCoordsFromDir() const{
    double l,b;

    //get the transformation matrix from celestial to galactic coordinates.
    Rotation celToGal(celToGal());
    
    //and do the transform to get the galactic celestial vector
    Hep3Vector pointingin(celToGal*m_dir);
    
    // pointingin is the galactic cartesian pointing vector,
    // we want to make this into l and b now.
    l = atan(pointingin.x()/pointingin.z());
    //b = atan(pointingin.y()/pointingin.z());
    b = asin(pointingin.y());
    
    l *= 360./M_2PI;
    b *= 360./M_2PI;
    
    //a serious kluge - this part needs further examination
    if(pointingin.z()<0){
        if(pointingin.x()>=0){
            l=180.+l;
        }else if(pointingin.x()<=0){
            l=-180.+l;
        }
    }
    return std::make_pair<double,double>(l,b);
}

std::pair<double,double> SkyDir::setCelCoordsFromDir() const{
    double ra,dec;

    //we now want to use the cartesian vector to get (ra, dec).
    ra = atan(m_dir.x()/m_dir.z());
    //m_dec = atan(m_dir.y()/m_dir.z());
    dec = asin(m_dir.y());
    
    ra *= 360./M_2PI;
    dec *= 360./M_2PI;
    return std::make_pair<double,double>(ra,dec);
}

double SkyDir::l ()const{
    std::pair<double,double> abc = setGalCoordsFromDir();
    return setGalCoordsFromDir().first;
}

double SkyDir::b ()const{
    return setGalCoordsFromDir().second;
}

double SkyDir::ra ()const{
    return setCelCoordsFromDir().first;
}

double SkyDir::dec ()const{
    return setCelCoordsFromDir().second;
}

Hep3Vector SkyDir::r ()const{
    return m_dir;
}

