#ifndef TAGGERHIT_HH
#define TAGGERHIT_HH


// Our own classes.
#include "AdfEvent/FadcDataWord.h"

namespace AncillaryData
{
  class TaggerHit
    {
    public:
		TaggerHit(){
			m_IsPedestalSubtracted = false;
			m_moduleId = 0;
			m_stripId = 0;
			m_layerId = 0;
			m_pulseHeight = 0;
		};

      TaggerHit(FadcDataWord fadcDw)
	{
	  m_IsPedestalSubtracted = false;
	  m_moduleId = fadcDw.getFadcModule();
	  m_stripId  = fadcDw.getFadcChannel();
	  m_layerId  = (m_stripId<384) ? 0 : 1 ;
	  m_pulseHeight=fadcDw.getFadcValue();
	}
      ~TaggerHit(){;}
      
      unsigned    getModuleId()     const {return m_moduleId;}
      void        setModuleId(unsigned moduleId) { m_moduleId = moduleId ; }
      unsigned    getStripId()      const {return m_stripId;}
      void        setStripId(unsigned stripId) { m_stripId = stripId; }
      unsigned    getPulseHeight()  const {return m_pulseHeight;}
      void setPulseHeight(unsigned pulseHeight) { m_pulseHeight = pulseHeight; }
      unsigned    getLayerId() const {return m_layerId;}
      void        setLayerId(unsigned layerId) { m_layerId = layerId; }
      bool        getPedestalSubtract()const {return m_IsPedestalSubtracted;}
      void        setPedestalSubtract(){m_IsPedestalSubtracted = true;}
      void        print(){std::cout<<"Tagger Module: "<<getModuleId()<<" ch: "<<getStripId()<<" ("<<getLayerId()<<") PH: "<<getPulseHeight()<<std::endl;}
    private:
      
      bool m_IsPedestalSubtracted;
      unsigned m_moduleId;
      unsigned m_layerId;
      unsigned m_stripId;
      unsigned m_pulseHeight;  
    };
}
#endif
  
