//  $Header$
//  Author: Toby Burnett
//  Project: graphics
//

#ifndef GRAPHICS_VECTOR_H
#define GRAPHICS_VECTOR_H


class SpacePoint;
class SpaceVector;
class Hep3Vector;
class HepPoint3D;
namespace gui {


class GraphicsVector {
// 
// Purpose: Define a special 3-vector, for use by the drawing programs only
// Go to this trouble since we want to have a float representation for
// economy of virtual memory, and to be in principle independent of other
public:
    GraphicsVector(float x=0, float y=0, float z=0):m_x(x), m_y(y), m_z(z){}


    // should define conversion constructors for all relevant vector or points    
    // conversion from a SpacePoint, for convenience since returned by transformations
#if 0 // D0
    GraphicsVector(const SpacePoint& v);

    GraphicsVector(const SpaceVector& v);
#endif
    GraphicsVector(const Hep3Vector & v);
    GraphicsVector(const HepPoint3D & v);

    float x()const{return m_x;}
    float y()const{return m_y;}
    float z()const{return m_z;}
    
    GraphicsVector& operator+=(const GraphicsVector&);
    GraphicsVector  operator-(const GraphicsVector&);
    GraphicsVector  operator+(const GraphicsVector&);
    double mag2()const; // 2-norm

private:
    float m_x, m_y, m_z;
};

inline GraphicsVector& GraphicsVector::operator += (const GraphicsVector & p) {
  m_x += p.x(); 
  m_y += p.y();
  m_z += p.z();
  return *this;
}
inline GraphicsVector GraphicsVector::operator - (const GraphicsVector & p) {
  return GraphicsVector(m_x-p.x(),m_y-p.y(),m_z-p.z());
}
inline GraphicsVector GraphicsVector::operator + (const GraphicsVector & p) {
  return GraphicsVector(m_x+p.x(),m_y+p.y(),m_z+p.z());
}
inline double GraphicsVector::mag2()const{return m_x*m_x+m_y*m_y+m_z*m_z;}

inline GraphicsVector operator * (float a, const GraphicsVector & p) {
  return GraphicsVector(a*p.x(), a*p.y(), a*p.z());
}

} //namespace gui
#endif //GRAPHICS_VECTOR_H
