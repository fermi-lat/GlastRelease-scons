//$Header$

/// Gaudi specific include files
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

// gui, display includes
#include "GuiSvc/IGuiSvc.h"
#include "gui/DisplayControl.h"
#include "gui/GuiMgr.h"
#include "CalRecon/CalDisplay.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "idents/VolumeIdentifier.h"
#include "CLHEP/Geometry/Transform3D.h"

using namespace Event;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
static const AlgFactory<CalDisplay>  Factory;
const IAlgFactory& CalDisplayFactory = Factory;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CalDisplay::CalDisplay(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator) 
{
    
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class CalRep : public gui::DisplayRep {
private:
	IDataProviderSvc* m_eventSvc;
	float m_xtalHeight;
    float m_calZtop;
    float m_calZbottom;

public:
    CalRep(IDataProviderSvc* eventSvc, float xtalHeight,
        float calZtop, float calZbottom)
		:m_eventSvc(eventSvc),m_xtalHeight(xtalHeight),
         m_calZtop(calZtop), m_calZbottom(calZbottom){}
    void update(){

		const Point p0(0.,0.,0.);
		CalXtalRecCol* cxrc = SmartDataPtr<CalXtalRecCol>(m_eventSvc,"/Event/CalRecon/CalXtalRecCol"); 
		if(cxrc){

// drawing red box for each log with a size proportional to energy deposition
			
			setColor("red");

			double emax = 0.;
			for (CalXtalRecCol::const_iterator it = cxrc->begin();
				it != cxrc->end(); it++){
				CalXtalRecData* recData = *it;

				double eneXtal = recData->getEnergy();
				if(eneXtal>emax)emax=eneXtal;
			}

			
			
			if(emax>0){

				for (CalXtalRecCol::const_iterator it = cxrc->begin();
						it != cxrc->end(); it++){
						CalXtalRecData* recData = *it;

						double eneXtal = recData->getEnergy();
						if(eneXtal>0.01*emax){


							Vector pXtal = recData->getPosition() - p0;
							int layer = (recData->getPackedId()).getLayer();
							double x = pXtal.x();
							double y = pXtal.y();
							double z = pXtal.z();
							double s = 0.45*m_xtalHeight*eneXtal/emax;
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

//		drawing the cross in the average position for each layer 
		CalClusterCol* cls = SmartDataPtr<CalClusterCol>(m_eventSvc,"/Event/CalRecon/CalClusterCol");
		if(cls){
			double s=0.1*m_xtalHeight;
			setColor("blue");
			CalCluster* cl = cls->getCluster(0);
			double energy_sum = cl->getEnergySum();
			const std::vector<double>& eneLayer = cl->getEneLayer();
			const std::vector<Vector>& posLayer = cl->getPosLayer();
			for( int l=0;l<8;l++){
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


//		drawing the center of the cluster		
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
            
            double dirX = (cl->getDirection()).x();
            double dirY = (cl->getDirection()).y();
            double dirZ = (cl->getDirection()).z();
        
            if(dirZ >= -1. && dirZ != 0.){

                double xTop = x+dirX*(m_calZtop-z)/dirZ;
                double yTop = y+dirY*(m_calZtop-z)/dirZ;
                double xBottom = x+dirX*(m_calZbottom-z)/dirZ;
                double yBottom = y+dirY*(m_calZbottom-z)/dirZ;

                moveTo(Point(xTop,yTop,m_calZtop));
                lineTo(Point(xBottom,yBottom,m_calZbottom));


            }

        }
	}
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
StatusCode CalDisplay::initialize()
{
    //Look for the gui service
    MsgStream log(msgSvc(), name());
    IGuiSvc* guiSvc = 0;
    StatusCode sc = service("GuiSvc", guiSvc);
    

    sc = service("GlastDetSvc", detSvc);

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
    if(!detSvc->getNumericConstByName(std::string("CALnLayer"), &value)) {
       log << MSG::ERROR << " constant " << " CALnLayer "<<" not defined" << endreq;
       return StatusCode::FAILURE;
    } else nLayers = value;
    if(!detSvc->getNumericConstByName(std::string("eLATTowers"), &value)) {
       log << MSG::ERROR << " constant " << " eLATTowers "<<" not defined" << endreq;
       return StatusCode::FAILURE;
    } else eLATTowers = value;
    if(!detSvc->getNumericConstByName(std::string("eTowerCAL"), &value)) {
       log << MSG::ERROR << " constant " << " eTowerCAL "<<" not defined" << endreq;
       return StatusCode::FAILURE;
    } else eTowerCAL = value;
    if(!detSvc->getNumericConstByName(std::string("eXtal"), &value)) {
       log << MSG::ERROR << " constant " << " eXtal "<<" not defined" << endreq;
       return StatusCode::FAILURE;
    } else eXtal = value;

    if(!detSvc->getNumericConstByName(std::string("CsIHeight"), &xtalHeight)) {
       log << MSG::ERROR << " constant " << " CsIHeight "<<" not defined" << endreq;
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
    
    
    
    //Ok, see if we can set up the display
    if (sc.isSuccess())  {
        //Set up the display rep for Clusters
        guiSvc->guiMgr()->display().add(
            new CalRep(eventSvc(),xtalHeight,calZtop,calZbottom), "Cal recon");
    }
    
    return sc;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
StatusCode CalDisplay::execute()
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    return sc;
}


