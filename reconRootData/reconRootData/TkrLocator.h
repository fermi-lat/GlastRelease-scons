#ifndef TKRLOCATOR_H
#define TKRLOCATOR_H

#include "TObject.h"
#include "TVector3.h"
#include "TVector2.h"

// ===================================
//          TkrLocator
// ===================================
//! Root object specifying the position & direction of a track
/*! Specifies the position and direction of a track at a certain locations 
    within the instrument */

class TkrLocator : public TObject          
{                                     

public : 
    //! define a list of locations...any more locations??
    typedef enum {
        VERTEX = 0,
        LASTTKRLAYER,
        CAL,
        ACD,
        OTHER
    } LocationType;

private :
  
  //! type of location (vertex, CAL, ...)
  LocationType m_type;  
  //! id of the TkrLocator
  UInt_t m_id;           
  
  //! x,y,z in mm
  TVector3* m_position;  
  //! px/pz, py/pz
  TVector2* m_direction;  

  // 2 symmetric matrices defined by 3 elements
  //! error in x
  Float_t m_sigma_x;        
  //! error in slope for x
  Float_t m_sigma_xslope;  
  // covariance for x
  Float_t m_cov_x;           

  //! error in y
  Float_t m_sigma_y;      
  //! error in slope for y
  Float_t m_sigma_yslope;  
  //! covariance for y
  Float_t m_cov_y;    


public : 

    // constructors
    // No TObjArray in this class - so we don't have to
    // worry about the default/non-default constructor
    // issue
    //! default constructor
    TkrLocator();
    //! constructor that allows one to specify the ID #
    TkrLocator(UInt_t id);

    //! destructor
    virtual ~TkrLocator();

    void Clean();
    void Create();

    //! access the ID # of this TkrLocator
    inline UInt_t getId() { return m_id; };
    //! set the ID of this TkrLocator
    inline void setId(UInt_t id) { m_id = id; };

    //! access the type, defined by the location 
    inline LocationType getType() { return m_type; };
    //! set the type
    inline void setType(LocationType loc) { m_type = loc; };

    //! access the position
    const TVector3& getPosition() { return *m_position; };
    //! set the position of this TkrLocator
    inline void setPosition(TVector3 *pos) { m_position = pos; };

    //! access the direction
    inline const TVector2& getSlope() { return *m_direction; };
    //! set the direction of this TkrLocator in X
    inline Double_t getSlopeX() { return m_direction->X(); };
    //! set the direction of this TkrLocator in Y
    inline Double_t getSlopeY() { return m_direction->Y(); };
    //! set the direction where we've provided both x and y
    inline void setDirection(TVector2 *dir) { m_direction = dir; };

    //! access the error in x
    inline Float_t getSigma_X() { return m_sigma_x; };
    //! set the error in x
    inline void setSigma_X(Float_t sigma) { m_sigma_x = sigma; };
    
    inline Float_t getSigma_Xslope() { return m_sigma_xslope; };
    inline void setSigma_Xslope(Float_t sigma) { m_sigma_xslope = sigma; };
    
    inline Float_t getCov_X() { return m_cov_x; };
    inline void setCov_X(Float_t cov) { m_cov_x = cov; };

    //! setup the 3 elements of the X matrix (x error, slope error, covariance)
    void setXmatrix(Float_t x, Float_t slope, Float_t cov);

    inline Float_t getSigma_Y() { return m_sigma_y; };
    inline void setSigma_Y(Float_t sigma) { m_sigma_y = sigma; };

    inline Float_t getSigma_Yslope() { return m_sigma_yslope; } ;
    inline void setSigma_Yslope(Float_t sigma) { m_sigma_yslope = sigma; };

    inline Float_t getCov_Y() { return m_cov_y; };
    inline void setCov_Y(Float_t cov) { m_cov_y = cov; };

    //! setup the 3 elements of the Y matrix (y error, slope error, covariance)
    void setYmatrix(Float_t y, Float_t slope, Float_t cov);

  ClassDef(TkrLocator,1) 
};
#endif
