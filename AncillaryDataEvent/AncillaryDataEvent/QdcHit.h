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
		  m_pulseHeight = 0.0;
	  } ;

      QdcHit(QdcDataWord qdcDw)
	{
	  m_IsPedestalSubtracted=false;
	  m_channel     = qdcDw.getQdcChannel();
	  m_module      = qdcDw.getQdcModule();
	  m_pulseHeight = static_cast<double>(qdcDw.getQdcValue());
	}
      ~QdcHit(){;}
      
      unsigned int    getQdcModule()       const {return m_module;}
      void            setQdcModule(unsigned int mo)   { m_module = mo; }

      unsigned int    getQdcChannel()      const {return m_channel;}
      void            setQdcChannel(unsigned int ch) { m_channel = ch; }

      double          getPulseHeight()      const {return m_pulseHeight;}
      void            setPulseHeight(double pulseHeight) { m_pulseHeight = pulseHeight; }

      void            SubtractPedestal(double pedestalValue, double rms)
	{
	  const double PHPS = getPulseHeight()-pedestalValue;
	  setPulseHeight(PHPS);
	  setPedestalSubtract();
	  setSigma(rms);
	}
      void setSigma(double rms){m_sigma = rms;}
      double getSigma() const {return m_sigma;}      

      bool       getPedestalSubtract() const {return m_IsPedestalSubtracted;}
      void       setPedestalSubtract(){m_IsPedestalSubtracted=true;}
      
      void        print(){std::cout<<" ------- QCD HIT Module: "<<getQdcModule()<<" Ch: "<<getQdcChannel()<<" PH "<<getPulseHeight()<<std::endl;}
    private:
      
      bool m_IsPedestalSubtracted;
      unsigned int m_module;
      unsigned int m_channel;
      double   m_pulseHeight;  
      double m_sigma;
    };
}
#endif
