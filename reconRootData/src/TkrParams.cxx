#include "reconRootData/TkrParams.h"
#include <iostream>

ClassImp(TkrParams)

TkrParams::TkrParams() 
{
    Clear();
}


TkrParams::TkrParams(Double_t ax, Double_t sx, Double_t ay, Double_t sy) {
    initialize(ax, sx, ay, sy);
}

TkrParams::TkrParams(const TkrParams& params)
{
    initialize(params);
}

void TkrParams::Clear(Option_t *option) {
    m_xpos = 0.0;
    m_xslope = 0.0;
    m_ypos = 0.0;
    m_yslope = 0.0;
}

void TkrParams::Print(Option_t *option) const {
    TObject::Print(option);
    using namespace std;

    cout << "(xpos, xslope) = (" << m_xpos << ", " << m_xslope << ")" << endl;
    cout << "(ypos, yslope) = (" << m_ypos << ", " << m_yslope << ")" << endl;

}

TkrParams& TkrParams::operator=(const TkrParams& params)
{
    initialize(params);

    return *this;
}

void TkrParams::initialize(Double_t ax, Double_t sx, Double_t ay, Double_t sy)
{
    m_xpos = ax;
    m_xslope = sx;
    m_ypos = ay;
    m_yslope = sy;
}  

void TkrParams::initialize(const TkrParams& params)
{
    m_xpos = params.getXPos();
    m_xslope = params.getXSlope();
    m_ypos = params.getYPos();
    m_yslope = params.getYSlope();
}  


const Double_t& TkrParams::operator() (const int &element) const
{
    return (*this)[element];
}

Double_t& TkrParams::operator() (const int &element)
{
    return (*this)[element];
}


const Double_t& TkrParams::operator[] (const int &element) const
{
    // Purpose and Method:  Element access using [] operator.
    // Output:  Element corresponding to index number

  if (element == 0 ) {
    return m_xpos;
  } else if (element == 1) {
    return m_xslope;
  } else if (element == 2) {
    return m_ypos;
  } else if (element == 3) {
    return m_yslope;
  } else {
    return m_xpos;
  }  
      
}

Double_t& TkrParams::operator[] (const int &element)
{
    if (element == 0 ) {
        return m_xpos;
    } else if (element == 1) {
        return m_xslope;
    } else if (element == 2) {
        return m_ypos;
    } else if (element == 3) {
        return m_yslope;
    } else {
        return m_xpos;
    }  
    
}

const TkrParams TkrParams::operator +(const TkrParams& A) const
{
    
    Double_t xpos = m_xpos + A.m_xpos;
    Double_t xslope = m_xslope + A.m_xslope;
    Double_t ypos = m_ypos + A.m_ypos;
    Double_t yslope = m_yslope + A.m_yslope;
    
    return TkrParams(xpos,xslope,ypos,yslope);
}

const TkrParams TkrParams::operator -(const TkrParams& A) const
{
    
    Double_t xpos = m_xpos - A.m_xpos;
    Double_t xslope = m_xslope - A.m_xslope;
    Double_t ypos = m_ypos - A.m_ypos;
    Double_t yslope = m_yslope - A.m_yslope;
    
    return TkrParams(xpos,xslope,ypos,yslope);
    
}

const Double_t TkrParams::operator *(const TkrParams& A) const
{
    
    return m_xpos * A.m_xpos + m_xslope * A.m_xslope
        + m_ypos * A.m_ypos + m_yslope * A.m_yslope;
    
}


