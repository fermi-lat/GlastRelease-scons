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
		  m_pulseHeight = 0;
	  } ;

      QdcHit(QdcDataWord qdcDw)
	{
	  m_IsPedestalSubtracted=false;
	  m_channel  = qdcDw.getQdcChannel();
	  m_pulseHeight=qdcDw.getQdcValue();
	}
      ~QdcHit(){;}
      
      unsigned    getQdcChannel()      const {return m_channel;}
      void        setQdcChannel(unsigned ch) { m_channel = ch; }
      unsigned    getPulseHeight()  const {return m_pulseHeight;}

      void setPulseHeight(unsigned pulseHeight) { m_pulseHeight = pulseHeight; }
      
      bool       getPedestalSubtract() const {return m_IsPedestalSubtracted;}
      void       setPedestalSubtract(){m_IsPedestalSubtracted=true;}
      void        print(){std::cout<<"QCD Ch: "<<getQdcChannel()<<" PH "<<getPulseHeight()<<std::endl;}
    private:
      
      bool m_IsPedestalSubtracted;
      unsigned m_channel;
      unsigned m_pulseHeight;  
    };
}
#endif
