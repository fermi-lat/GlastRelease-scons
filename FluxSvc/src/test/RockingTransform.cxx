
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ToolFactory.h"

#include "FluxSvc/IRegisterSource.h"
#include "FluxSvc/ISpectrumFactory.h"
#include "FluxSvc/IFluxSvc.h"

//needed?
//#include "geometry/Point.h"
//#include "geometry/Vector.h"
//#include "geometry/CoordTransform.h"

#include "CLHEP/Vector/Rotation.h"

/** @class RockingTransform
 *  @brief Register and implement the rocking angle tool
 *  
 *   @author Sean Robinson
 *   $Header$
 */
class RockingTransform : public AlgTool, virtual public IRegisterSource {
public:
    
    RockingTransform( const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~RockingTransform() { }
    
    /// implement to define interface: will be called from FluxSvc
    StatusCode registerMe(IFluxSvc* );

    /// return the rotation from the specified rocking angles.
    Rotation rockingRotation(double time);

    ///micro-utility to display the rotation matrix on the screen.
    void displayRotation(Rotation rot);
    
};


// Static factory for instantiation of algtool objects
static ToolFactory<RockingTransform> s_rockingFactory;
const IToolFactory& RockingTransformFactory = s_rockingFactory;

// Standard Constructor
RockingTransform::RockingTransform(const std::string& type, 
                         const std::string& name, 
                         const IInterface* parent)
                         : AlgTool( type, name, parent ) {
    
    // Declare additional interface
    declareInterface<IRegisterSource>(this);
}


StatusCode RockingTransform::registerMe(IFluxSvc* fsvc) 
{
    MsgStream  log(msgSvc(), name());
    log << MSG::INFO << "Rocking Transform Tool Registered, testing..." << endreq;
    //static RemoteSpectrumFactory<MapSpectrum> factory(fsvc);
    displayRotation(rockingRotation(0.));
    return StatusCode::SUCCESS;
} 

Rotation RockingTransform::rockingRotation(double time){
    //Purpose:  return the rotation to correct for satellite rocking.
    //Input:  Current time
    //Output:  3x3 rocking-angle transformation matrix.
    Rotation gal;   
    //and here we construct the rotation matrix
    double zenithPhase = 0.0;//should be input?  or time-based?;
    double offZenith = 0.0;//input?  time-based?
    gal.rotateZ(zenithPhase).rotateX(offZenith);
    
    return gal;
}


void RockingTransform::displayRotation(Rotation rot){
    
    std::cout << "{" << rot.xx() << ' ' << rot.xy() << ' ' << rot.xz() << std::endl
        << rot.yx() << ' ' << rot.yy() << ' ' << rot.yz() << std::endl
        << rot.zx() << ' ' << rot.zy() << ' ' << rot.zz() << "}" << std::endl;
}