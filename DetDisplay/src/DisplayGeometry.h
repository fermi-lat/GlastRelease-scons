/** @file DisplayGeometry.h 
     @brief Declatation of class  DisplayGeometry
  $Header$
*/

#ifndef DisplayGeometry_h
#define DisplayGeometry_h

#include <vector>
#include <map>
#include <fstream>
#include <strstream>

#include "geometry/CoordTransform.h"
#include "GlastSvc/GlastDetSvc/IGeometry.h"
#include "idents/VolumeIdentifier.h"
class Shape;
class Medium;
/**
This class instantiates the Glast  geometry for display, from detModel
*/
class DisplayGeometry : public IGeometry
{
public:
    /**
    Create the object. 
    Expect the function members shape, push, pop, and id to be called in proper sequence.
  */
    DisplayGeometry();
    ~DisplayGeometry();
    /**

        @param s type of the shape
        @param id vector of unsigned ints (maybe null)
        @param name
        @param material
        @param params vector with the six transformation parameters, followed by 3 or so dimensions
        @return tell caller whether to skip subvolumes or not
    */
    virtual VisitorRet pushShape(ShapeType s, const UintVector& id, 
                                 std::string name, std::string material, 
                                 const DoubleVector& params, VolumeType type,
                                 SenseType sense);
    
    //* called to signal end of nesting */
    virtual void popShape();

    void printStats(std::ostream& out);

    /// current ID
    idents::VolumeIdentifier getId()const; 

    typedef std::map<std::string, std::pair<int,double> > MaterialSummary;
    const MaterialSummary & materials()const{ return m_matSum;}

private:
    typedef std::vector<std::string> StringVector;

    /// create a shape, box or whatever. If material is "", it is an enclosing volume
    /// parameters (x,y,z for a Box) are the dimensions of the shape.
    /// the type specifies how it is to be built:
    ///   Composite: will be filled with sub-volumes
    ///   XStack: interior volumes will be stacked in the local x-direction
    ///   Simple or Sensitive: leaf on the tree, not containing further vols
    bool shape(ShapeType s, std::string name, std::string material, const DoubleVector& params, VolumeType type, SenseType sense);
    
    /// push a transformation onto the stack
    void push(double x, double y, double z, double rx, double ry, double rz);
    
    
    /// push a new id onto the idstack (beware, may be more than one per level)
    virtual void id( std::string name, double value);

    //! keep track of the current transformation
    std::vector<CoordTransform> m_Tstack;

    //! stack of boxes that must be hierarchical
    std::vector<Shape* > m_vols;

    //! stack of material names
    std::vector<std::string> m_materialStack;

    //! stack of CompositeMedium objects
    std::vector<Medium*> m_mediumStack;
    
    //! vector of ids to describe current object
    UintVector m_idValues;

    /// should be equivalent
    idents::VolumeIdentifier m_volId;

    //! vector of the number of ids for this geometry level: 0,1, or even 2
    UintVector m_idcount;

    // this accounts for total volume associated with each material: index by material name,
    // save count and total volume

    MaterialSummary m_matSum;
    /// choice mode for traversing geometry
    std::string m_mode;

};
#endif
