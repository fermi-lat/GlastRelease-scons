#ifndef POSITION_H
#define POSITION_H
#include "detModel/Sections/Boundable.h"
#include "detModel/Sections/Units.h"
#include "detModel/Sections/IdField.h"
#include "detModel/Utilities/Purge.h"
#include <vector>

namespace detModel{

  class SectionsVisitor;
  class Units;
  class Volume;
  /**
   * This class represents a generic positioned volume
   * @author R.Giannitrapani
   */
  class Position : public Boundable{
  public:
  
    /// A standard constructor
    Position():volume(0),volumeRef(""){units = new Units;}

    /// In the destructor we need to purge the idFields vector and the units
    virtual ~Position(){purge(idFields); delete units;}
    virtual void Accept(SectionsVisitor*) = 0;
    virtual void AcceptNotRec(SectionsVisitor*) =0;

    /**
     * This methods sets and gets the volume and volumeref
     */ 
    void setVolumeRef(std::string pname){volumeRef=pname;}
    void setVolume(Volume* pvolume){volume=pvolume;}
    std::string getVolumeRef()const{return volumeRef;}
    Volume* getVolume()const{return volume;}

    /**
     *  This methods add an idField and retrive the vector of idField
     */ 
    void addIdField(IdField* pidField){idFields.push_back(pidField);}
    std::vector <IdField*> getIdFields(){return idFields;}

  private:
    /** @link aggregation */
    std::vector <IdField*> idFields;
    Units* units;
    /** @link aggregation
     *   @supplierCardinality 1 */
    Volume * volume;
    std::string volumeRef;
  
  };
}
#endif //POSITION_H
