// GraphicsVector.cpp

#include "gui/GraphicsVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Geometry/Point3D.h"

using namespace gui;

// implement conversion constructors here

GraphicsVector::GraphicsVector(const Hep3Vector& v)
    : m_x(v.x()), m_y(v.y()), m_z(v.z()) {}

GraphicsVector::GraphicsVector(const HepPoint3D& v)
    : m_x(v.x()), m_y(v.y()), m_z(v.z()) {}
