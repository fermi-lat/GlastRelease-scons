//$Header$

/// Gaudi specific include files
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

// gui, display includes
#include "GuiSvc/IGuiSvc.h"
#include "gui/GuiMgr.h"
#include "CalDisplay.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "idents/VolumeIdentifier.h"
#include "CLHEP/Geometry/Transform3D.h"

using namespace Event;

static const AlgFactory<CalDisplay>  Factory;
const IAlgFactory& CalDisplayFactory = Factory;

CalDisplay::CalDisplay(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator) {}


void CalRep::update()

// Purpose and method: this function extracts the reconstructed calorimeter
//                     data from TDS classes and draw the graphic
//                     representation using methods of gui::DisplayRep class
//
// TDS Inputs: CalXtalRecCol and CalClustersCol 


{
    
    const Point p0(0.,0.,0.);

    
    // get the pointer to CalXtalRecCol object in TDS
    CalXtalRecCol* cxrc = SmartDataPtr<CalXtalRecCol>(m_eventSvc,
        EventModel::CalRecon::CalXtalRecCol); 

    // if pointer is not zero - draw the reconstructed xtal data
    if(cxrc){
        
        // drawing red box for each log with a size
        // proportional to energy deposition
        
        setColor("red");
        
        double emax = 0.; // reset maximum energy per crystal
        
        // to find maximum energy per crystal
        for (CalXtalRecCol::const_iterator it = cxrc->begin();
        it != cxrc->end(); it++){

            // get poiner to the reconstructed data for individual crystal
            CalXtalRecData* recData = *it;
            
            // get reconstructed energy in the crystal
            double eneXtal = recData->getEnergy();

            // if energy is bigger than current maximum - update the maximum 
            if(eneXtal>emax)emax=eneXtal;
        }
        
        
        // if maximum crystal energy isn't zero - start drawing 
        if(emax>0){
            

            // loop over all crystals in reconstructed collection
            // to draw red boxes
            for (CalXtalRecCol::const_iterator it = cxrc->begin();
            it != cxrc->end(); it++){
                
                // get poiner to the reconstructed data for individual crystal
                CalXtalRecData* recData = *it;
                
                // get reconstructed energy in the crystal
                double eneXtal = recData->getEnergy();

                
                // draw only crystals containing more than 1% of maximum energy
                if(eneXtal>0.01*emax){
                    
                    // get the vector of reconstructed position
                    Vector pXtal = recData->getPosition() - p0;
                    
                    // get reconstructed coordinates
                    double x = pXtal.x();
                    double y = pXtal.y();
                    double z = pXtal.z();

                    
                    // calculate the half size of the box, 
                    // taking the 90% of crystal half height
                    // as the size corresponding to the maximum energy
                    double s = 0.45*m_xtalHeight*eneXtal/emax;

                    // drawing the box with the center at the reconstructed
                    // position (x,y,x)
                    moveTo(Point(x-s, y-s, z-s));
                    lineTo(Point(x+s, y-s, z-s));
                    lineTo(Point(x+s, y-s, z+s));
                    lineTo(Point(x-s, y-s, z+s));
                    lineTo(Point(x-s, y-s, z-s));
                    moveTo(Point(x-s, y+s, z-s));
                    lineTo(Point(x+s, y+s, z-s));
                    lineTo(Point(x+s, y+s, z+s));
                    lineTo(Point(x-s, y+s, z+s));
                    lineTo(Point(x-s, y+s, z-s));
                    moveTo(Point(x-s, y-s, z-s));
                    lineTo(Point(x-s, y+s, z-s));
                    moveTo(Point(x+s, y-s, z-s));
                    lineTo(Point(x+s, y+s, z-s));
                    moveTo(Point(x-s, y-s, z+s));
                    lineTo(Point(x-s, y+s, z+s));
                    moveTo(Point(x+s, y-s, z+s));
                    lineTo(Point(x+s, y+s, z+s));
                }
            }
        }
    }
    
    // drawing the cross in the average position for each layer 

    //  get pointer to the cluster reconstructed collection
    CalClusterCol* cls = SmartDataPtr<CalClusterCol>(m_eventSvc,
        EventModel::CalRecon::CalClusterCol);

    
    // if pointer is not zero, start drawing
    if(cls){

        // set the cross half size to 10% of crystal height
        double s=0.1*m_xtalHeight;
        setColor("blue"); // set cross color to blue

        // get pointer to the cluster 0 - the only one exiting now
        CalCluster* cl = cls->getCluster(0); 

        // get total energy in the calorimeter
        double energy_sum = cl->getEnergySum();

        // get vector of layer energies
        const std::vector<double>& eneLayer = cl->getEneLayer();

        // get layer positions
        const std::vector<Vector>& posLayer = cl->getPosLayer();

        // draw only if there is some energy in the calorimeter        
        if(energy_sum > 0){        

            // loop over calorimeter layers
            for( int l=0;l<8;l++){

                // if energy in this layer is not zero - draw blue cross at
                // the average reconstructed position for this layer
                if (eneLayer[l]>0){
                    double x=(posLayer[l]).x();
                    double y=(posLayer[l]).y();
                    double z=(posLayer[l]).z();
                    moveTo(Point(x-s, y, z));
                    lineTo(Point(x+s, y, z));
                    moveTo(Point(x, y-s, z));
                    lineTo(Point(x, y+s, z));
                    moveTo(Point(x, y, z-s));
                    lineTo(Point(x, y, z+s));
                    
                }
            }
            
            
            // drawing the center of the cluster		
            setColor("green");
            double x = (cl->getPosition()).x();
            double y = (cl->getPosition()).y();
            double z = (cl->getPosition()).z();
            moveTo(Point(x-s, y, z));
            lineTo(Point(x, y, z+s));
            lineTo(Point(x+s, y, z));
            lineTo(Point(x, y, z-s));
            lineTo(Point(x-s, y, z));
            lineTo(Point(x, y+s, z));
            lineTo(Point(x+s, y, z));
            lineTo(Point(x, y-s, z));
            lineTo(Point(x-s, y, z));
            moveTo(Point(x, y-s, z));
            lineTo(Point(x, y, z+s));
            lineTo(Point(x, y+s, z));
            lineTo(Point(x, y, z-s));
            lineTo(Point(x, y-s, z));		
            

            // drawing the reconstructed shower direction
            // as a green line
            double dirX = (cl->getDirection()).x();
            double dirY = (cl->getDirection()).y();
            double dirZ = (cl->getDirection()).z();
            
            // non display for non-physical or horizontal direction
            if(dirZ >= -1. && dirZ != 0.){
                
                
                // calculate x and y coordinates for the beginning and the end
                // of line in the top and bottom calorimeter layers
                double xTop = x+dirX*(m_calZtop-z)/dirZ;
                double yTop = y+dirY*(m_calZtop-z)/dirZ;
                double xBottom = x+dirX*(m_calZbottom-z)/dirZ;
                double yBottom = y+dirY*(m_calZbottom-z)/dirZ;
                
                //draw shower direction
                moveTo(Point(xTop,yTop,m_calZtop));
                lineTo(Point(xBottom,yBottom,m_calZbottom));
                
                
            }
        }        
    }
}



StatusCode CalDisplay::initialize()
{
    //Look for the gui service
    MsgStream log(msgSvc(), name());
    
    
    StatusCode sc = service("GlastDetSvc", detSvc);
    
    if(sc.isFailure())
    {
        log << MSG::ERROR << "GlastDetSvc could not be found" <<endreq;
        return sc;
    }
    double xtalHeight; 
    int nLayers;
    int eLATTowers;
    int eTowerCAL;
    int eXtal;
    
    double value;
    if(!detSvc->getNumericConstByName(std::string("CALnLayer"), &value))
    {
        log << MSG::ERROR << " constant " << " CALnLayer "
            <<" not defined" << endreq;
        return StatusCode::FAILURE;
    } else nLayers = int(value);
    if(!detSvc->getNumericConstByName(std::string("eLATTowers"), &value))
    {
        log << MSG::ERROR << " constant " << " eLATTowers "
            <<" not defined" << endreq;
        return StatusCode::FAILURE;
    } else eLATTowers = int(value);
    if(!detSvc->getNumericConstByName(std::string("eTowerCAL"), &value))
    {
        log << MSG::ERROR << " constant " << " eTowerCAL "
            <<" not defined" << endreq;
        return StatusCode::FAILURE;
    } else eTowerCAL = int(value);
    if(!detSvc->getNumericConstByName(std::string("eXtal"), &value))
    {
        log << MSG::ERROR << " constant " << " eXtal "
            <<" not defined" << endreq;
        return StatusCode::FAILURE;
    } else eXtal = int(value);
    
    if(!detSvc->getNumericConstByName(std::string("CsIHeight"),&xtalHeight)) 
    {
        log << MSG::ERROR << " constant " << " CsIHeight "
            <<" not defined" << endreq;
        return StatusCode::FAILURE;
    } 
    
    
    int layer=0;
    idents::VolumeIdentifier topLayerId;
    topLayerId.append(eLATTowers);
    topLayerId.append(0);
    topLayerId.append(0);
    topLayerId.append(eTowerCAL);
    topLayerId.append(layer);
    topLayerId.append(layer%2);
    topLayerId.append(0);
    topLayerId.append(eXtal);
    topLayerId.append(0);
    
    HepTransform3D transfTop;
    detSvc->getTransform3DByID(topLayerId,&transfTop);
    Vector vecTop = transfTop.getTranslation();
    
    layer=nLayers-1;
    idents::VolumeIdentifier bottomLayerId;
    bottomLayerId.append(eLATTowers);
    bottomLayerId.append(0);
    bottomLayerId.append(0);
    bottomLayerId.append(eTowerCAL);
    bottomLayerId.append(layer);
    bottomLayerId.append(layer%2);
    bottomLayerId.append(0);
    bottomLayerId.append(eXtal);
    bottomLayerId.append(0);
    
    HepTransform3D transfBottom;
    detSvc->getTransform3DByID(bottomLayerId,&transfBottom);
    Vector vecBottom = transfBottom.getTranslation();
    
    
    float calZtop = vecTop.z();
    float calZbottom = vecBottom.z();
    
    
    IGuiSvc* guiSvc = 0;
    sc = service("GuiSvc", guiSvc);
    
    //Ok, see if we can set up the display
    if (sc.isSuccess())  {
        //Set up the display rep for Clusters
        guiSvc->guiMgr()->display().add(
            new CalRep(eventSvc(),xtalHeight,calZtop,calZbottom), "Cal recon");
    }
    
    return sc;
}



