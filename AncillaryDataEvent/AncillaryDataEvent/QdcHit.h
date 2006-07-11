#ifndef QDCHIT_HH
#define QDCHIT_HH


// Our own classes.
#include "AdfEvent/QdcDataWord.h"

namespace AncillaryData
{
  class QdcHit
    {
    public:
      
	  QdcHit() {
		  m_IsPedestalSubtracted = false;
		  m_channel = 0;
		  m_module = 0;
		  m_pulseHeight = 0;
	  } ;

      QdcHit(QdcDataWord qdcDw)
	{
	  m_IsPedestalSubtracted=false;
	  m_channel  = qdcDw.getQdcChannel();
	  m_module  = qdcDw.getQdcModule();
	  m_pulseHeight=qdcDw.getQdcValue();
	}
      ~QdcHit(){;}
      
      unsigned int    getQdcModule()       const {return m_module;}
      void            setQdcModule(unsigned int mo)   { m_module = mo; }

      unsigned int    getQdcChannel()      const {return m_channel;}
      void            setQdcChannel(unsigned int ch) { m_channel = ch; }

      unsigned int    getPulseHeight()      const {return m_pulseHeight;}
      void            setPulseHeight(unsigned int pulseHeight) { m_pulseHeight = pulseHeight; }
      
      bool       getPedestalSubtract() const {return m_IsPedestalSubtracted;}
      void       setPedestalSubtract(){m_IsPedestalSubtracted=true;}
      
      void        print(){std::cout<<"QCD HIT Module: "<<getQdcModule()<<" Ch: "<<getQdcChannel()<<" PH "<<getPulseHeight()<<std::endl;}
    private:
      
      bool m_IsPedestalSubtracted;
      unsigned int m_module;
      unsigned int m_channel;
      unsigned int m_pulseHeight;  
    };
}
#endif
