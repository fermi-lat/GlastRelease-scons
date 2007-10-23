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
    enum axisDir{xDir, yDir, zDir};
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
    virtual void visitTrap(Trap*);
    virtual void visitSphere(Sphere*);
    virtual void visitPosXYZ(PosXYZ*);
    virtual void visitAxisMPos(AxisMPos*);
    virtual void visitIdField(IdField*);
    virtual void visitSeg(Seg*);
    virtual void insertVolume(Volume*);

    /// This method returns the full map of positioned volume pointers
    const std::map<idents::VolumeIdentifier, const PositionedVolume*>* 
    getVolMap() const;

    /// Retrieve the HepTransform3D from an ID.  ID should be global
    /// (using all fields, from world volume on down) and transform
    /// will also be relative to world volume
    /// Returns false if fails to find the id
    bool getTransform3DByID(idents::VolumeIdentifier, HepGeom::Transform3D*);

    /// Retrieve the HepTransform3D from an ID.  ID should be relative
    /// to top volume.  However transform will be relative to world volume. 
    /// See getTopTransformById to retrieve transform relative to top volume
    /// Returns false if fails to find the id
    bool getTransform3DByTopID(idents::VolumeIdentifier, HepGeom::Transform3D*);

    /// Retrieve the HepTransform3D from an ID.  ID should be relative
    /// to top volume.  Transform will also be relative to top volume.
    /// Returns false if fails to find the id
    bool getTopTransform3DByTopID(idents::VolumeIdentifier, HepGeom::Transform3D*);

    /// Retrieve the HepTransform3D from an ID.  ID should be complete,
    /// but transform will be relative to top volume.
    /// Returns false if fails to find the id
    bool getTopTransform3DByID(idents::VolumeIdentifier, HepGeom::Transform3D*);

    /// Retrieve the shape parameters and type from an ID
    /// Returns false if fails to find the id
    bool getShapeByID(idents::VolumeIdentifier, 
                      std::string* s, std::vector<double>* params);

    /// Retrieve the shape parameters and type from an ID where ID
    /// is relative to top volume ID.
    /// Returns false if fails to find the id
    bool getShapeByTopID(idents::VolumeIdentifier, 
                      std::string* s, std::vector<double>* params);

    /// Return ID of top volume relative to world
    idents::VolumeIdentifier getIDPrefix() const {return m_IDPrefix;}

    /// Return transform of top volume relative to world
    const HepGeom::Transform3D& getTransform3DPrefix() const 
    {return m_transformPrefix;}

    /** 
	This method returns a PositionedVolume pointer from an ID (a null 
        pointer	if the ID does not exist.  ID must be full ID (all fields
        starting with world volume)
    */
    const PositionedVolume* getPositionedVolumeByID(idents::VolumeIdentifier) 
      const;
    typedef std::map<idents::VolumeIdentifier, const PositionedVolume*> PVmap;

    /**
      Given a sample identifier and specification of fields which may be
      wildcarded,  append identifiers which match, ordered by coordinate
      (X, Y, or Z) and sign (increasing or decreasing) specified, to
      caller-supplied vector.  This routine does not clear the vector
      first.
     */
    unsigned orderSensitive(std::vector<idents::VolumeIdentifier>& vols,
                            const idents::VolumeIdentifier& sample,
                            const std::vector<bool> wild,
                            axisDir coord, 
                            bool increasing = true) const;

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
    CLHEP::Hep3Vector m_actualPos;
    /// This is the actual rotation of the volume
    CLHEP::HepRotation m_actualRot;
    /// This is the map of PositionedVolume pointers indexed by ids
    PVmap m_volMap;

    /// a secondary vector to preserve the order
    IdVector m_idvec;

    /// Set to true when our top volume is encountered
    bool m_topSeen;

    /// Prefix needed to convert from identifier for a volume relative
    /// to topVolume and identifier for the same volume relative to
    /// the full geometry (using top volume specified by section element)
    idents::VolumeIdentifier m_IDPrefix;

    HepGeom::Transform3D m_transformPrefix;
    HepGeom::Transform3D m_inverseTransformPrefix; // avoid recomputing

    Gdd*  m_gdd;
  };

}
#endif //IDMAPBUILDER







