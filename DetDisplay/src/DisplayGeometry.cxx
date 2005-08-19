/** @file DisplayGeometry.cxx
     @brief  implementation of class DisplayGeometry
  $Header$
*/
#include "DisplayGeometry.h"

#include "idents/VolumeIdentifier.h"

// geometry only
#include "geometry/Box.h"
#include "geometry/CoordTransform.h"
#include "geometry/Vector.h"
#include "geometry/Tube.h"

// gismo 
#include "World.h"

// for the display
#include "gui/GuiMgr.h"
#include "gui/DisplayControl.h"
#include "gui/DisplayRep.h"
#include "geomrep/TubeRep.h"
#include "geomrep/BoxRep.h"
#include <typeinfo>
#include <iomanip>

#include <sstream>
#include <cassert>
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  special Medium subclasses for overriding graphics, detector setup
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class GlastCompositeMedium : public CompositeMedium {
public:
    GlastCompositeMedium(Medium* mom, Shape* vol, const char * mat)
        :CompositeMedium(mom, vol, mat){}
    
    void createDetectorView(gui::DisplayRep& v){
        v.setColor("white");
        Medium::createDetectorView(v);
        for(iterator it=begin(); it !=end(); ++it )
            (*it)->createDetectorView(v.nested());
        // this winds up hiding all but top level box
        v.hide(); 
    }
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class GlastMedium : public Medium {
public:
    GlastMedium(Medium* mom, Shape* vol, const char * mat, bool sensitive)
        :Medium(mom, vol, mat),m_sensitive(sensitive){}
    
    void createDetectorView(gui::DisplayRep& v) {

        const std::type_info& t = typeid(volume());
        v.setColor(m_sensitive? "cyan" : "grey" );
        if( t==typeid(Box) )    v.append(BoxRep(volume() ));
        else if (t==typeid(Tube) )v.append(TubeRep(volume() ));
        else {
        }
    }
    bool m_sensitive;
};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DisplayGeometry::DisplayGeometry()
{
    
    // prime the stack with identity for convenience
    m_Tstack.push_back(CoordTransform());
    m_materialStack.push_back(std::string("Vacuum"));

    m_mediumStack.push_back( World::instance());
    
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DisplayGeometry::~DisplayGeometry()
{
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
IGeometry::VisitorRet 
  DisplayGeometry::pushShape(ShapeType s, const UintVector& idvec, 
                           std::string name, std::string material, 
                           const DoubleVector& params, 
                           VolumeType type, SenseType sense)
{
    push(params[0],params[1],params[2],params[3],params[4],params[5]);
    IGeometry::DoubleVector::const_iterator v=params.begin();
    IGeometry::DoubleVector shapepars;
    v+=6;
    std::copy(v,params.end(),std::back_inserter<IGeometry::DoubleVector>(shapepars));
    for( UintVector::const_iterator u=idvec.begin(); u!=idvec.end(); ++u) id("",*u); 
    bool abort = shape(s,name,material,shapepars,type, sense);
    return abort ? AbortSubtree :  More;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayGeometry::push(double x, double y,double z, double rx, double ry, double rz)
{
    Rotation R;
    // note that this is the global to local transformation: signs are reversed.
    if( rx!=0) R.rotateX(rx*M_PI/180.);
    if( ry!=0) R.rotateY(ry*M_PI/180.);
    if( rz!=0) R.rotateZ(rz*M_PI/180.);
     m_Tstack.push_back(CoordTransform(R,Vector(x,y,z))*m_Tstack.back());
    m_idcount.push_back(0);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool DisplayGeometry::shape(ShapeType shape_type, std::string name, std::string material, 
                          const DoubleVector& params, VolumeType volType,
                            SenseType sense)
{
    double volume=0;
    bool abort = false;
    ::Shape* s=0;
    if( shape_type==Box) {
        double x=params[0], y=params[1], z=params[2];
        volume=x*y*z;
        s= new ::Box(x,y,z); 
    
    }else if(shape_type==Tube) {
        double dz=params[0], rmin=params[1], rmax=params[2];
        s = new ::Tube(dz, rmin, rmax);
        volume= M_PI*rmax*rmax*dz;
    }
    assert(s);


    // check to see if mother vol has material to be reduced.
    std::string momMatName = m_materialStack.back();
    if( !momMatName.empty() && momMatName != std::string("Vacuum")) {
        m_matSum[momMatName].second-= volume;
    }
    // translate material name
    std::string glastMaterial = material;
    if(0 ==&material || material==std::string("Vacuum") ) glastMaterial = "vacuum";
    
    // special to abort further traversal of detailed tree for Si and CsI
    if(  name.find("SiLadders")!= std::string::npos)  { glastMaterial = "Si"; abort = true; }
    else if( name == "CsIDetector"){ 
        glastMaterial = "CsI"; abort=true;  }

    // add to material summary
    std::pair<int,double> entry= m_matSum[glastMaterial];
    ++entry.first;
    entry.second += volume;
    m_matSum[glastMaterial] = entry;

    
    // put into Medium hierarchy
    Medium* newmed;
    if( volType !=  Simple && ! abort)
        newmed = new GlastCompositeMedium(m_mediumStack.back(), s, glastMaterial.c_str());
    else newmed = new GlastMedium(m_mediumStack.back(), s, glastMaterial.c_str(),
        sense==posSensitive || sense==intSensitive || abort);
    
    m_mediumStack.push_back(newmed);
    m_vols.push_back(s);
    m_materialStack.push_back(glastMaterial);
    std::stringstream full_name;
    if(name.find("SiLayerBox")!= std::string::npos || 
        name.find("SiLadders")!= std::string::npos ) {
        // special to put in names expected by KalParticle
        unsigned int top = m_idValues.back();
        full_name<< ( top==1? "top strip": "bot strip");
    } else {
        
        full_name << name; 
        
        for( UintVector::iterator ui=m_idValues.begin(); ui !=m_idValues.end(); ++ui) 
            full_name << "/"<<+*ui;
    }
    full_name <<'\0';

    newmed->setTitle(full_name.str().c_str());
    if( sense==posSensitive || sense==intSensitive || name=="oneTower") {
        //TODO: take care of CsI hit display?
    }

    // now transform the medium and detector if any
    newmed->transform(m_Tstack.back());
    return abort;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayGeometry::printStats(std::ostream& out)
{
    int total = 0;
    out << " Summary of created volumes " << std::endl;
    double vol=0; 
    out << std::setw(15) << "material"  << std::setw(8)  << "count"  << std::setw(12) << "vol(cm^3)" << std::endl;
    out << std::setw(15) << "--------"  << std::setw(8)  << "-----"  << std::setw(12) << "---------" << std::endl;
    for( MaterialSummary::const_iterator it = m_matSum.begin(); it != m_matSum.end(); ++it) {
        if( (*it).first.empty() ) continue;
        out << std::setw(15) << (*it).first 
            << std::setw(8)  << (*it).second.first 
            << std::setw(12) <<  static_cast<int>((*it).second.second/1e3) << std::endl;
        total += (*it).second.first ;
        vol   += (*it).second.second;
    }
    out << std::setw(23)  << "------"  << std::setw(12) << "------" <<std::endl;
    out << std::setw(15) << "Totals" << std::setw(8) << total << std::setw(12)<< static_cast<int>(vol/1e3) 
        << std::endl << std::endl;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayGeometry::popShape()
{
    m_Tstack.pop_back();
    if( !m_vols.empty() ) m_vols.pop_back();
    m_mediumStack.pop_back();
    m_materialStack.pop_back();
    
    unsigned int count = m_idcount.back(); m_idcount.pop_back();
    while(count--){
        m_idValues.pop_back();
    }
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayGeometry::id( std::string name, double value)
{
    unsigned int this_id = static_cast<unsigned int>(value);
    m_idValues.push_back(this_id);
    //m_idNames.push_back(name);
    ++ m_idcount.back(); 
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
idents::VolumeIdentifier DisplayGeometry::getId()const
{
    idents::VolumeIdentifier id;
    for(UintVector::const_iterator i=m_idValues.begin(); i!=m_idValues.end(); ++i) id.append(*i);
    return id;
}
