#ifndef IDMAPBUILDER
#define IDMAPBUILDER
#include <vector>
#include <map>

#include "detModel/Management/SectionsVisitor.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

namespace detModel{
  class PositionedVolume;

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
    IDmapBuilder(std::string);     
      
    virtual ~IDmapBuilder();
  
    /**
     * This is the visitor for the GDDsectionsContainer */
    virtual void visitGdd(Gdd*);
  
    /**
     * This is the visitor for the Section 
     */
    virtual void visitSection(Section*);

    /**
     * This is the visitor for the Ensemble 
     */
    virtual void visitEnsemble(Ensemble*);

    /**
     * This is the visitor for the Box 
     */
    virtual void visitBox(Box*);

    /**
     * This is the visitor for the Tube 
     */
    virtual void visitTube(Tube*);

    /**
     * This is the visitor for the PosXYZ 
     */
    virtual void visitPosXYZ(PosXYZ*);

    /**
     * This is the visitor for the AxisMPos 
     */
    virtual void visitAxisMPos(AxisMPos*);

    /**
     * This is the visitor for the IdField 
     */
    virtual void visitIdField(IdField*);

    /**
     * This is the visitor for the Seg 
     */
    virtual void visitSeg(Seg*);


    virtual void insertVolume(Volume*);

    std::string m_actualVolume;
    std::string m_actualID;
    
    Hep3Vector m_actualPos;
    HepRotation m_actualRot;

    
    std::map<std::string, PositionedVolume*> m_volMap;
  };

}
#endif //IDMAPBUILDER







