// $Header$

// Include files
#include <iostream>
#include "geometry/Vector.h"
#include "Event/Utilities/SkyDir.h"
#include "geometry/CoordTransform.h"

SkyDir::SkyDir(double param1, double param2, coordSystem inputType){
    if(inputType == GALACTIC){
        m_l = param1;
        m_b = param2;
        
        if (m_l==0.){m_l+=0.000000000001;}  //to fix divide-by-zero errors
        
        //here we construct the cartesian galactic vector
        Vector gamgal(sin(m_l*M_2PI/360.)*cos(m_b*M_2PI/360.) , sin(m_b * M_2PI/360.) , cos(m_l*M_2PI/360.)*cos(m_b*M_2PI/360.));
        
        //get the transformation matrix from galactic to celestial coordinates.
        Rotation galToCel=celToGal().inverse();
        
        //and do the transform to get the cartesian celestial vector
        m_dir = galToCel*gamgal;
        
        setCelCoordsFromDir();
        
    }else if(inputType == CELESTIAL){
        m_ra = param1;
        m_dec = param2;
        
        if (m_ra==0.){m_ra+=0.000000000001;}  //to fix divide-by-zero errors
        
        //here we construct the cartesian celestial vector
        m_dir = Vector(sin(m_ra*M_2PI/360.)*cos(m_dec*M_2PI/360.) , sin(m_dec*M_2PI/360.) , cos(m_ra*M_2PI/360.)*cos(m_dec*M_2PI/360.));        
        
        setGalCoordsFromDir();
        
    }else{
        //improper coordinate system declaration - default things and say so.
        std::cout << "Improper coordinate System declaration in SkyDir" << std::endl;
        m_l = 0;
        m_b = 0;
        m_ra = 0;
        m_dec = 0;
        
        m_dir = Vector(0,0,0);
    }
}

SkyDir::SkyDir(Vector dir):
m_dir(dir){
    setCelCoordsFromDir();
    setGalCoordsFromDir();
}

Rotation SkyDir::celToGal(){
    //gal is the matrix which rotates (cartesian)celestial coordiantes into (cartesian)galactic ones
    Rotation gal;
    double degsPerRad = 180./M_PI;
    gal.rotateZ(-282.25/degsPerRad).rotateX(-62.6/degsPerRad).rotateZ(33./degsPerRad);
    return gal; 
}


void SkyDir::setGalCoordsFromDir(){
    
    //get the transformation matrix from celestial to galactic coordinates.
    Rotation celToGal(celToGal());
    
    //and do the transform to get the galactic celestial vector
    Vector pointingin(celToGal*m_dir);
    
    // pointingin is the galactic cartesian pointing vector,
    // we want to make this into l and b now.
    m_l = atan(pointingin.x()/pointingin.z());
    //b = atan(pointingin.y()/pointingin.z());
    m_b = asin(pointingin.y());
    
    m_l *= 360./M_2PI;
    m_b *= 360./M_2PI;
    
    //a serious kluge - this part needs further examination
    if(pointingin.z()<0){
        if(pointingin.x()>=0){
            m_l=180.+m_l;
        }else if(pointingin.x()<=0){
            m_l=-180.+m_l;
        }
    }
}

void SkyDir::setCelCoordsFromDir(){
    //we now want to use the cartesian vector to get (ra, dec).
    m_ra = atan(m_dir.x()/m_dir.z());
    //m_dec = atan(m_dir.y()/m_dir.z());
    m_dec = asin(m_dir.y());
    
    m_ra *= 360./M_2PI;
    m_dec *= 360./M_2PI;
    
    //a serious kluge - this part needs further examination
    /*if(m_dir.z()<0){
        if(m_dir.x()>=0){
            m_ra=180.+m_ra;
        }else if(m_dir.x()<=0){
            m_ra=-180.+m_ra;
        }
    } */     
}

double SkyDir::l(){
    return m_l;
}

double SkyDir::b(){
    return m_b;
}

double SkyDir::ra(){
    return m_ra;
}

double SkyDir::dec(){
    return m_dec;
}

Vector SkyDir::r(){
    return m_dir;
}

