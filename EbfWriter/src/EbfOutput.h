#ifndef EBF_OUTPUT_H
#define EBF_OUTPUT_H


#include <stdio.h>

class EbfAcdData;
class EbfTkrData;
class EbfCalData;
class EbfGemData;









class EbfOutput
{
  public:

//   inline EbfOutput ()  { return; }
   EbfOutput::EbfOutput () : LdfFormat(false) { return; }
 
   void         setLdfFormat(bool flg) { LdfFormat=flg; return;}      
   int          open   (const char    *fileName,
                        unsigned int maxEvtSize);

   unsigned int format (const EbfAcdData *acd,
                        const EbfCalData *cal,
                        const EbfTkrData *tkr,
                        const EbfGemData *gem,
                        unsigned int *mcInfo);
   
   void         writeMC (unsigned int *mcInfo, int size);
   void         print  ()  const;
   char * write  (bool writeEbf, unsigned int &length, unsigned int *TdsBuffer);
   char * getData(unsigned int &dataSize) {dataSize=m_curEvtSize; return (char *)m_evtBuffer;}
   unsigned int flush  ();
   unsigned int close  ();
   inline int   setPrint  (int printFlag) 
   {
       int old = m_print;
       m_print = printFlag;
       return old;
   }
           
   
   inline void setNumEvtsOut(unsigned int n) { m_numEvtsOut = n; }
   inline unsigned int numEvtsOut () const { return m_numEvtsOut; }
   inline unsigned int numEvtsIn  () const { return m_numEvtsIn;  }   
   
   
   ~EbfOutput ();
   
  /*
   | The unit of all size parameters is 32 - bit words
  */ 
  private:
    int                  m_print; /*!< Print flag                           */ 
    FILE                   *m_fp; /*!< The file handle                      */
    unsigned int    m_maxEvtSize; /*!< The maximum size of one event        */
    unsigned int    m_curEvtSize; /*!< The current size in the event buffer */
    unsigned int    m_totEvtSize; /*!< The total size of the events written */
    unsigned int     m_numEvtsIn; /*!< The number of  input events          */
    unsigned int    m_numEvtsOut; /*!< The number of output events          */
    unsigned int    *m_evtBuffer; /*!< The event buffer                     */
    unsigned int   m_circBuffOff; /*!< Offset in circular buffer           */
    unsigned int    *m_evtHead;   /*!< Point to event header               */
    unsigned int    *m_evtDescriptor; /*!< Event Descriptor   */
    unsigned int    *m_evtEnd;     /*!< Pointer to end of event */
    bool            LdfFormat;     /*! < Flag for LDF Format */
};




#endif




