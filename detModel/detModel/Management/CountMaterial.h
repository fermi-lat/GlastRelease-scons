#ifndef COUNTMATERIAL_H
#define COUNTMATERIAL_H
#include "detModel/Management/SectionsVisitor.h"
#include "detModel/Sections/Box.h"
#include "detModel/Sections/Tube.h"
#include <map>

namespace detModel{

  /**
   *  This visitor counts the number of volumes that has
   *  a particular material assigned and fill a map
   *  @author R.Giannitrapani
   */
  class CountMaterial : public SectionsVisitor {

  public:
    /// This is the constructor
    CountMaterial(std::string);
    
    /**
     * This is the visitor for the Gdd 
     */
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
    virtual void visitIdField(IdField*){};
    /**
     * This is the visitor for the Seg
     */
    virtual void visitSeg(Seg*){};

    /// This method return the number of volumes given a material name
    double getNumber(std::string name){return m_matNumber[name];}

    
    /// This method dumps to std output the matNumber table
    void printTable();

    /** 
       This is the map of double indicized by the name of material:
       it is double and not int since in future it will contain the volume
       in mm3 of material.
    */
    std::map<std::string, double> m_matNumber;

    /// The name of the topVolume
    std::string m_actualVolume;    
  };

}
#endif //COUNTMATERIAL_H







