#ifndef SCALERHIT_HH
#define SCALERHIT_HH


// Our own classes.
#include "AdfEvent/ScalerDataWord.h"

namespace AncillaryData
{
  class ScalerHit
    {
    public:
      
      ScalerHit() {
	m_channel = 0;
	m_scalerValue = 0;
      } ;
      
      ScalerHit(ScalerDataWord scalerDw, unsigned int channel)
	{
	  m_channel  = channel;
	  m_scalerValue=scalerDw.getScalerValue();
	}
      ~ScalerHit(){;}
      
      unsigned int    getScalerChannel()      const {return m_channel;}
      void            setScalerChannel(unsigned int ch) { m_channel = ch; }

      unsigned int    getScalerValue()      const {return m_scalerValue;}
      void            setScalerValue(unsigned int scalerValue) { m_scalerValue = scalerValue;}
      
      void        print(){std::cout<<"Scaler HIT  Ch: "<<getScalerChannel()<<" Value "<<getScalerValue()<<std::endl;}
    private:
      

      unsigned int m_channel;
      unsigned int m_scalerValue;
    };
}
#endif
