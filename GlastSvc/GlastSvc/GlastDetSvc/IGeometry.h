
#ifndef GLASTDET_IGEOMETRY_H
#define GLASTDET_IGEOMETRY_H
#include <vector>
#include <string>

/** @class IGeometry
* @brief Abstract class that defines callbacks from a traversal of the 
* detector geometry tree.
*
* $Header$
*/
class IGeometry {
public:

  typedef std::vector<double> DoubleVector;
  typedef std::vector<unsigned int>UintVector;
  enum VolumeType{ Simple, posSensitive, intSensitive, Composite, Xstack, Ystack, Zstack };
  enum ShapeType{ Box, Tube };
  enum VisitorRet { More, AbortSubtree};

  /** Notification of a new shape

      @param s type of the shape
      @param id vector of unsigned ints (maybe null)
      @param name of the shape
      @param material  material name, possibly meaningful in detModel
      @param params vector with the six transformation parameters, followed by 3 or so dimensions
      @return indicate whether traversal should continue normally or skip 
              current subtree.  Would be nice to also support total abort
        */
  virtual VisitorRet pushShape(ShapeType s, const UintVector& id, 
                               std::string name, std::string material, 
                               const DoubleVector& params, 
                               VolumeType type)=0;
    
  //* called to signal end of nesting */
  virtual void popShape()=0;

  /** Default implementation of getMode() will cause GeometryVisitor to
      use mode it gets from the visitorMode property.  */
  virtual std::string getMode() {return std::string();}
  virtual ~IGeometry(){}

protected:
    IGeometry(){};
};



#endif
