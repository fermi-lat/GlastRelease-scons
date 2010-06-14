#include "reconRootData/TkrLocator.h"

ClassImp(TkrLocator)

TkrLocator::TkrLocator() {
    // default constructor
  m_type = VERTEX;       // type of location ( vertex, ....)
  m_id = 0;              // id of the object (needed ??)
  
  m_position = 0;        // position in mm
  m_direction = 0;          // direction of the track at m_position

  m_sigma_x = 0.0f;        // error in x
  m_sigma_xslope = 0.0f;   // error in slope for x
  m_cov_x = 0.0f;            // covariance for x

  m_sigma_y = 0.0f;        // error in y
  m_sigma_yslope = 0.0f;   // error in slope for y
  m_cov_y = 0.0f;            // covariance for y
}

TkrLocator::TkrLocator(UInt_t id) {
    // non-default constructor
    // must be used or our pointers won't be initialized correctly
  m_id = id;
  m_type = VERTEX;            // type of location ( vertex, ....)
  
  m_sigma_x = 0.0f;        // error in x
  m_sigma_xslope = 0.0f;   // error in slope for x
  m_cov_x = 0.0f;            // covariance for x

  m_sigma_y = 0.0f;        // error in y
  m_sigma_yslope = 0.0f;   // error in slope for y
  m_cov_y = 0.0f;            // covariance for y

  m_position = 0;
  m_direction = 0;
}

TkrLocator::~TkrLocator() {
    Clean();
}

void TkrLocator::Clean() {
    if (m_position) {
        delete m_position;
        m_position = 0;
    }
    
    if (m_direction) {
        delete m_direction;
        m_direction = 0;
    }

    m_sigma_x = 0.0f;
    m_sigma_xslope = 0.0f;
    m_cov_x = 0.0f;

    m_sigma_y = 0.0f;
    m_sigma_yslope = 0.0f;
    m_cov_y = 0.0f;
}

void TkrLocator::Create() {
}

void TkrLocator::setXmatrix(Float_t x, Float_t slope, Float_t cov) {
    m_sigma_x = x;            // error in x
    m_sigma_xslope = slope;   // error in slope for x
    m_cov_x = cov;            // covariance for x
}

void TkrLocator::setYmatrix(Float_t y, Float_t slope, Float_t cov) {
    m_sigma_y = y;        // error in y
    m_sigma_yslope = slope;   // error in slope for y
    m_cov_y = cov;          // covariance for y
}
