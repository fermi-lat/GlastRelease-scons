#include "Event/Recon/CalRecon/GcrReconClasses.h"
#include <string>

//-----------------------------------------------------------------------------------------------------------------
void Event::GcrXtal::initialize(idents::CalXtalId xtalId, double pathLength, double closestFaceDist, int crossedFaces, Point entryPoint, Point exitPoint)

{
    m_xtalId   = xtalId;
    m_pathLength = pathLength;
    m_closestFaceDist = closestFaceDist ;
    m_crossedFaces  = crossedFaces;
    m_entryPoint = entryPoint ;
    m_exitPoint  = exitPoint;
}

//-----------------------------------------------------------------------------------------------------------------
void Event::GcrXtal::writeOut(MsgStream& log) const
{
    log << MSG::DEBUG;
    if (log.isActive() ) 
    {
        int it=m_xtalId.getTower();
        int il=m_xtalId.getLayer();
        int ic=m_xtalId.getColumn();
        log << "---> writeOut GcrXtal tow / lay / col = " << it << " " << il << " " << ic << endreq;
        //log << "---> x / y /z = " << m_xtalData->getPosition().x() << " " << m_xtalData->getPosition().y() << " " << m_xtalData->getPosition().z() << endreq;
        log << "---> pathLength = " << m_pathLength << endreq;
    }
}
//-----------------------------------------------------------------------------------------------------------------

void Event::GcrXtal::getReadableXedFaces(int xedFaces, int& inFace, int& outFace){

   bool m_debugging = false;
   
   if(m_debugging)
       std::cout<< "GcrXtal BEGIN getReadableXedFaces(), xedFAces = " << xedFaces<< std::endl;
   
   typedef enum {
        XFACE_ZTOP,
        XFACE_ZBOT,
        XFACE_XLEFT,
        XFACE_XRIGHT,
        XFACE_YLEFT,
        XFACE_YRIGHT
    } XFACE_BITPOS;

   int n, subXedFaces;
   // int m;
   std::string faceIn, faceOut;
   
   float mGuess;
  
   for(n=0; n<=5 ; n++){
     subXedFaces = int(xedFaces - std::pow(2.,double(n)));
     
     if(m_debugging)
         std::cout<< "subXedFaces= " << subXedFaces << std::endl;

     mGuess = std::log(double(subXedFaces))/std::log(double(2));
     
     if(m_debugging)
        std::cout<< "mGuess= " << mGuess << std::endl;

     if ((mGuess-int(mGuess)) == 0.){
       
       if(m_debugging)
           std::cout<< "m found, m= " << mGuess << std::endl;

       inFace=n;
       outFace=int(mGuess);

       
       if(m_debugging)
           std::cout<< "inFace,outFace found= " << inFace<< ","<< outFace << std::endl;
      
      
       break;  
     }
       
   
   }// end of for(n=0..5)
   
    if(m_debugging){
          
        if(inFace == XFACE_ZTOP){
            std::cout<< "XFACE_ZTOP"<< std::endl;
            faceIn =  "XFACE_ZTOP";
        } else if(inFace == XFACE_ZBOT){
            std::cout<< "XFACE_ZBOT"<< std::endl;
            faceIn =  "XFACE_ZBOT";
        }else if(inFace == XFACE_XLEFT){
            std::cout<< "XFACE_XLEFT"<< std::endl;
            faceIn =  "XFACE_XLEFT";
        }else if(inFace == XFACE_XRIGHT){
            std::cout<< "XFACE_XRIGHT"<< std::endl;
            faceIn =  "XFACE_XRIGHT";
        }else if(inFace == XFACE_YLEFT){
            std::cout<< "XFACE_YLEFT"<< std::endl;
            faceIn =  "XFACE_YLEFT";
        }else if(inFace == XFACE_YRIGHT){
            std::cout<< "XFACE_YRIGHT"<< std::endl;
            faceIn =  "XFACE_YRIGHT";
        }else{
            std::cout<< "DEFAULT"<< std::endl;
            faceIn =  "DEFAULT";
        }

        if(outFace == XFACE_ZTOP){
            std::cout<< "XFACE_ZTOP"<< std::endl;
            faceOut =  "XFACE_ZTOP";
        } else if(outFace == XFACE_ZBOT){
            std::cout<< "XFACE_ZBOT"<< std::endl;
            faceOut =  "XFACE_ZBOT";
        }else if(outFace == XFACE_XLEFT){
            std::cout<< "XFACE_XLEFT"<< std::endl;
            faceOut =  "XFACE_XLEFT";
        }else if(outFace == XFACE_XRIGHT){
            std::cout<< "XFACE_XRIGHT"<< std::endl;
            faceOut =  "XFACE_XRIGHT";
        }else if(outFace == XFACE_YLEFT){
            std::cout<< "XFACE_YLEFT"<< std::endl;
            faceOut =  "XFACE_YLEFT";
        }else if(outFace == XFACE_YRIGHT){
            std::cout<< "XFACE_YRIGHT"<< std::endl;
            faceOut =  "XFACE_YRIGHT";
        }else{
            std::cout<< "DEFAULT"<< std::endl;
            faceOut =  "DEFAULT";

        }
        
        std::cout<<"faceIn:"<< faceIn << std::endl;
        std::cout<<"faceOut:"<< faceOut << std::endl;

      }//end of if (m_debugging)


       //std::cout<< "faceIn,faceOut found= " << faceIn<< ","<< faceOut << std::endl;

       
     
   if(m_debugging)
       std::cout<< "GcrXtal END getReadableXedFaces()" << std::endl;


}




//-----------------------------------------------------------------------------------------------------------------
 std::ostream& Event::GcrXtal::fillStream( std::ostream& s ) const 
{ 
    s << "pathLength ="   << m_pathLength    << "\n";
  return s; 
}


//-----------------------------------------------------------------------------------------------------------------
void Event::GcrTrack::initialize(Vector direction, Vector dirError, Point calEntryPoint)
{
    m_direction   = direction;
    m_dirError = dirError;
    m_calEntryPoint = calEntryPoint;
}

