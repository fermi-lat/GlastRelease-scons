#ifndef GCRSelectClasses_H
#define GCRSelectClasses_H

#include "geometry/Point.h"
#include "geometry/Vector.h"

//#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/GcrReconClasses.h"

#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/MsgStream.h"

//static const CLID& CLID_CalMipTrackVecCol = InterfaceID("CalMipTrackVecCol", 1, 0);

/**   
* @class GcrSelectedXtal
*
* 
*/

//-----------------------------------------------------------------------------------------------------------------
namespace Event 
{
    //-----------------------------------------------------------------------------------------------------------------
    // Define a GcrXtal class which will be used in GCRSelectAlg
    class GcrSelectedXtal: public Event::GcrXtal, virtual public ContainedObject 
    {
        private:
            //Event::CalXtalRecData* m_xtalData;
            double                 m_rawEnergy;
            double                   m_corrEnergy;
            int                   m_selectGrade;
            

        public:
            GcrSelectedXtal(){};
              
            GcrSelectedXtal(idents::CalXtalId xtalId, double rawEnergy, double pathLength, double corrEnergy, int selectGrade, double closestFaceDist, int crossedFaces, Point entryPoint, Point exitPoint) : 
                        GcrXtal(xtalId, pathLength, closestFaceDist, crossedFaces, entryPoint, exitPoint), m_rawEnergy(rawEnergy), m_corrEnergy(corrEnergy), m_selectGrade(selectGrade){};

            virtual ~GcrSelectedXtal() {};

            void                   initialize(idents::CalXtalId xtalId, float rawEnergy, float pathLength, float corrEnergy, int selectGrade, double closestFaceDist, int crossedFaces, Point entryPoint, Point exitPoint);

            void                   setRawEnergy (double rawEnergy) {m_rawEnergy = rawEnergy ;}
            void                   setCorrEnergy  (double corrEnergy)         {m_corrEnergy  = corrEnergy;}
            void                   setSelectGrade  (int selectGrade)         {m_selectGrade  = selectGrade;}
           
            double                 getRawEnergy () const                               {return m_rawEnergy    ;}
            double                 getCorrEnergy () const                                {return m_corrEnergy     ;}
            int                   getSelectGrade () const                          {return m_selectGrade    ;}

          
            /// Utilities 
            void writeOut(MsgStream& log) const; 
            std::ostream& fillStream( std::ostream& s ) const;
    };
    

    // Define a vector of GcrSelectedXtals
    typedef std::vector<GcrSelectedXtal> GcrSelectedXtalsVec;

    // Define a vector of GcrSelectedXtals
    typedef ObjectVector<GcrSelectedXtal>      GcrSelectedXtalsCol;
    

    
    
    class GcrSelectVals: public DataObject
    {
    
            private:
        
            int m_inferedZ;
            int m_acdZ;
            int m_interactionParams;
            unsigned int m_gcrOBFStatusWord;



        public:
            GcrSelectVals(){};

            GcrSelectVals(int inferedZ, int acdZ, int interactionParams, unsigned int gcrOBFStatusWord) : 
                        m_inferedZ(inferedZ), m_acdZ(acdZ), m_interactionParams(interactionParams), m_gcrOBFStatusWord(gcrOBFStatusWord) {};

            ~GcrSelectVals() {};

            void                   initialize(int inferedZ, int acdZ, int interactionParams, unsigned int gcrOBFStatusWord);

            void                   setInferedZ (int inferedZ)         {m_inferedZ  = inferedZ;}
            void                   setAcdZ  (int acdZ)         {m_acdZ  = acdZ;}
            void                   setInteractionParams  (int interactionParams)         {m_interactionParams  = interactionParams;}
            void                   setGcrOBFStatusWord (unsigned int gcrOBFStatusWord){m_gcrOBFStatusWord  = gcrOBFStatusWord;}
           
            int                    getInferedZ ()  const                        {return m_inferedZ    ;}
            int                    getAcdZ ()      const                    {return m_acdZ    ;}
            double                 getInteractionParams ()  const                        {return m_interactionParams    ;}
            unsigned int           getGcrOBFStatusWord ()  const  {return m_gcrOBFStatusWord;}

           
            /// Utilities 
            //void writeOut(MsgStream& log) const; 
            //std::ostream& fillStream( std::ostream& s ) const;
    
    
    
    };

    
    
    
}
#endif
