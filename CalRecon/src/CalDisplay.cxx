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
#include "CalRecon/CsIClusters.h"
#include "CalRecon/CalRecLogs.h"
#include "CalRecon/CalGeometrySvc.h"
#include "CalRecon/CalDisplay.h"


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
	CalRecLogs** m_pp_crl;
	CsIClusterList** m_pp_cls;
	float m_logheight;
    float m_calZtop;
    float m_calZbottom;

public:
    CalRep(CalRecLogs** pp_crl, CsIClusterList** pp_cls, float logheight,
        float calZtop, float calZbottom)
		:m_pp_crl(pp_crl),m_pp_cls(pp_cls),m_logheight(logheight),
         m_calZtop(calZtop), m_calZbottom(calZbottom){}
    void update(){

		const Point p0(0.,0.,0.);
		CalRecLogs* crl = *m_pp_crl;
		if(crl){

// drawing red box for each log with a size proportional to energy deposition
			
			setColor("red");

			int nLogs = crl->num();
			double emax = 0.;
			for (int jlog = 0; jlog < nLogs ; jlog++) {
				CalRecLog* recLog = crl->Log(jlog);
				double eneLog = recLog->energy();
				if(eneLog>emax)emax=eneLog;
			}
			if(emax>0){
				for (jlog = 0; jlog < nLogs ; jlog++) {
					CalRecLog* recLog = crl->Log(jlog);
					double eneLog = recLog->energy();
					if(eneLog>0.01*emax){
						Vector pLog = recLog->position() - p0;
						double x = pLog.x();
						double y = pLog.y();
						double z = pLog.z();
						double s = 0.45*m_logheight*eneLog/emax;
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
		CsIClusterList* cls = *m_pp_cls;
		if(cls){
			double s=0.1*m_logheight;
			setColor("blue");
			CsICluster* cl = cls->Cluster(0);
			double energy_sum = cl->energySum();
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
			double x = (cl->position()).x();
			double y = (cl->position()).y();
			double z = (cl->position()).z();
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
            
            double dirX = (cl->direction()).x();
            double dirY = (cl->direction()).y();
            double dirZ = (cl->direction()).z();
        
            if(dirZ != 0.){

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
    IGuiSvc* guiSvc = 0;
    StatusCode sc = service("GuiSvc", guiSvc);

    m_crl = 0;
    m_cls = 0;
	
     sc = service("CalGeometrySvc", m_CalGeo);
	 float logheight = m_CalGeo->logHeight();
     float layerheight = m_CalGeo->layerHeight();
     int nlayers = m_CalGeo->numLayers()*m_CalGeo->numViews();
     float Z0 = m_CalGeo->Z0();
     float calZtop = Z0+(nlayers-1)*layerheight;
     float calZbottom = Z0;



    //Ok, see if we can set up the display
    if (sc.isSuccess())  {
	//Set up the display rep for Clusters
	guiSvc->guiMgr()->display().add(
        new CalRep(&m_crl,&m_cls,logheight,calZtop,calZbottom), "Cal recon");
    }
    
    return sc;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
StatusCode CalDisplay::execute()
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());


	m_crl = SmartDataPtr<CalRecLogs>(eventSvc(),"/Event/CalRecon/CalRecLogs"); 
	if (m_crl == 0){ sc = StatusCode::FAILURE;
	
	        log << MSG::ERROR << "CalDisplay failed to access CalRecLogs" << endreq;
			return sc;
    } else { m_crl->setCalDisplay(this);}
    
	
	m_cls  = SmartDataPtr<CsIClusterList>(eventSvc(),"/Event/CalRecon/CsIClusterList");

	if (m_cls == 0){ sc = StatusCode::FAILURE;
	
	        log << MSG::ERROR << "CalDisplay failed to access CsIClusterList" << endreq;
			return sc;
    } else {m_cls->setCalDisplay(this);}
    
    return sc;
}


