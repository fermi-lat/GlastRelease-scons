#ifndef IDMAPBUILDER
#define IDMAPBUILDER
#include <vector>
#include <map>

#include "detModel/Management/SectionsVisitor.h"
#include "idents/VolumeIdentifier.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Geometry/Transform3D.h"

namespace detModel{
  class PositionedVolume;
  class Volume;

  /**
   *  This visitor build a map of PositionedVolume indicized by ID.
   *  The method insertVolume(Volume*) contains the rules to add
   *  a particular volume to this map. For now it adds all the sensitive
   *  volumes (so shapes). In the future this class will become abstract
   *  and a concrete subclass will need to implement this insertVolume
   *  method. Also the map is public to access it easily, but it will 
   *  be changed to private with usual get method
   *  @author R.Giannitrapani
   */
  class IDmapBuilder : public SectionsVisitor {
    
  public:
    /**
     * The constructor has a string parameter; it is the name of the
     * root volume in GDD to represent along with all the volumes
     * contained in it. If it is the empty string the root volume is the
     * topVolume of the section.  */
    IDmapBuilder(std::string="");     
      
    virtual ~IDmapBuilder();
  
    virtual void visitGdd(Gdd*);
    virtual void visitSection(Section*);
    virtual void visitEnsemble(Ensemble*);
    virtual void visitBox(Box*);
    virtual void visitTube(Tube*);
    virtual void visitPosXYZ(PosXYZ*);
    virtual void visitAxisMPos(AxisMPos*);
    virtual void visitIdField(IdField*);
    virtual void visitSeg(Seg*);
    virtual void insertVolume(Volume*);

    /// This method return the full map of positioned volume pointers
    const std::map<idents::VolumeIdentifier, const PositionedVolume*>* getVolMap() const;

    /// Retrive the HepTransform3D from an ID 
    bool getTransform3DByID(idents::VolumeIdentifier, HepTransform3D*);


    /** 
	This method return a PositionedVolume pointer from an ID (a null pointer if
	the ID does not exist
    */
    const PositionedVolume* getPositionedVolumeByID(idents::VolumeIdentifier) const ;
    typedef std::map<idents::VolumeIdentifier, const PositionedVolume*> PVmap;

    //! make a summary on the stream
    void summary(std::ostream& );
    //! access to the container
    PVmap::const_iterator begin()const;
    PVmap::const_iterator end()const;
    size_t size()const;

    //! return pointer to a vector of ids, in order created
    typedef std::vector<idents::VolumeIdentifier> IdVector;
    const IdVector* getIdVector()const { return &m_idvec;}

  private:
    /// This string contains the topVolume used for the visit
    std::string m_actualVolume;
    /// This is the identifiers used during the m_volMap building
    idents::VolumeIdentifier m_actualID;    
    /// This is the actual position of the volume
    Hep3Vector m_actualPos;
    /// This is the actual rotation of the volume
    HepRotation m_actualRot;
    /// This is the map of PositionedVolume pointers indicized by ids
    PVmap m_volMap;

    /// a secondary vector to preserve the order
    IdVector m_idvec;

    Gdd*  m_gdd;
  };

}
#endif //IDMAPBUILDER







