#ifndef POSITIONEDVOLUME_H
#define POSITIONEDVOLUME_H

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"


namespace detModel{
  class Volume;

  /**
     This class represent a positioned volume; it contains a pointer 
     to a Volume (geometric and logical info), together with a rotation
     and translation (it uses CLHEP classes).

     @author R.Giannitrapani
   */
  class PositionedVolume {
    
  public:
    
    PositionedVolume():m_volume(0){};

    /// Set and get methods for the Volume pointer
    void setVolume(Volume* vol){m_volume = vol;};
    const Volume* getVolume(){return m_volume;};

    /// Set and get methods for the translation vector
    void setTranslation(Hep3Vector t){m_translation = t;};
    Hep3Vector getTranslation(){return m_translation;};

    /** 
       Set and get methods for the rotation matrix (to be changed 
       in order to return a reference or pointer       
    */
    void setRotation(HepRotation r){m_rotation = r;};
    HepRotation getRotation(){return m_rotation;};

  private:
    
    /// This is the Volume pointer of this positioned volume
    Volume* m_volume;
    /// The 3D translation
    Hep3Vector m_translation;
    /// The 3D rotation
    HepRotation m_rotation;
  };
}

#endif
