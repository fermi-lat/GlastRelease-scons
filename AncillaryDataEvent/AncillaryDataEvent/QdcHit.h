#ifndef QDCHIT_HH
#define QDCHIT_HH


// Our own classes.
#include "AdfEvent/QdcDataWord.h"

namespace AncillaryData
{
  class QdcHit
    {
    public:
      
      QdcHit(QdcDataWord qdcDw)
	{
	  m_IsPedestalSubtracted=0;
	  m_channel  = qdcDw.getQdcChannel();
	  m_pulseHeight=qdcDw.getQdcValue();
	}
      ~QdcHit(){;}
      
      unsigned    getQdcChannel()      const {return m_channel;}
      unsigned    getPulseHeight()  const {return m_pulseHeight;}
      
      bool       getPedestalSubtract(){return m_IsPedestalSubtracted;}
      void       setPedestalSubtract(){m_IsPedestalSubtracted=true;}
      void        print(){std::cout<<"QCD Ch: "<<getQdcChannel()<<" PH "<<getPulseHeight()<<std::endl;}
    private:
      
      bool m_IsPedestalSubtracted;
      unsigned m_channel;
      unsigned m_pulseHeight;  
    };
};
#endif
