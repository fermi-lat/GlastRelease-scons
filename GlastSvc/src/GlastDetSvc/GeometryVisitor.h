#ifndef GEOMETRYVISITOR_H
#define GEOMETRYVISITOR_H

#include "detModel/Management/SectionsVisitor.h"
#include <vector>

#include "GlastSvc/GlastDetSvc/IGeometry.h"
namespace detModel {class IdField; }

/**
* This is a concrete implementation of a "non recursive" sectionsVisitor that traverses
*  the geometry.
* It does no analysis, passing stuff on to a helper class interface IGeometry.
*/
class GeometryVisitor : public detModel::SectionsVisitor {
    
public:
    
    GeometryVisitor( IGeometry& geom ); 
    
    virtual void visitEnsemble(detModel::Ensemble*);
    
    virtual void visitBox(detModel::Box*);

    virtual void visitTube(detModel::Tube*);

    virtual void visitPosXYZ(detModel::PosXYZ*);
    
    virtual void visitAxisMPos(detModel::AxisMPos*);
    
    // following must be declared, but are not used
    virtual void visitGdd(detModel::Gdd*){}
    virtual void visitSection(detModel::Section*){};
    virtual void visitIdField(detModel::IdField*){}
    virtual void visitSeg(detModel::Seg*){}
    
private:
    //! typedefs for convenience
    typedef std::vector<detModel::IdField*> IdVector;
    typedef std::vector<float> FloatVector;

    void cacheId( detModel::IdField* id, int step); // note: should be const pointer!
    
    //! the object to call back, passed in with the constructor
    IGeometry& m_geom;
    IGeometry::DoubleVector m_params;

    void cacheTransform(double x=0, double y=0, double z=0, double rx=0, double ry=0, double rz=0);
    void pushTransform();

    IGeometry::UintVector m_idvec;
    void pushId();
};

#endif //GISMOGEOMETRYVISITOR_H
