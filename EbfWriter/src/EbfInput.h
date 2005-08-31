#ifndef EBF_INPUT_H
#define EBF_INPUT_H


#include <stdio.h>

class EbfAcdData;
class EbfTkrData;
class EbfCalData;
class EbfGemData;


class EbfInput
{
  public:

   inline EbfInput () { return; }
 
         
   int          open   (const char    *fileName,
                        unsigned int maxEvtSize,
                        bool ldfFormat);

   unsigned int format (const EbfAcdData *acd,
                        const EbfCalData *cal,
                        const EbfTkrData *tkr,
                        const EbfGemData *gem);

   void         print  ()  const;
   unsigned int read  (bool ldfFormat);
   unsigned int flush  ();
   unsigned int close  ();
   void         parse  (EbfTkrData *tkr, 
                        EbfCalData *cal,
                        EbfGemData *gem, 
                        EbfAcdData *acd,
                        EbfCalConstants *calCon);
   inline int   setPrint  (int printFlag) 
   {
       int old = m_print;
       m_print = printFlag;
       return old;
   }
           
   
   inline unsigned int numEvtsOut () const { return m_numEvtsOut; }
   inline unsigned int numEvtsIn  () const { return m_numEvtsIn;  }   
   
   
   ~EbfInput ();
   
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
    unsigned int    *m_evtHeader; /*!< Pointer to Event Header              */
    unsigned int     m_nReadEvt;
};




#endif




