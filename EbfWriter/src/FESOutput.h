#ifndef FES_OUTPUT_H
#define FES_OUTPUT_H


#include <stdio.h>


class EbfAcdData;
class EbfTkrData;
class EbfCalData;
class EbfGemData;

#define FES_FILEHEAD_SIZE 98
#define FES_RTYP_HEADER 0
#define FES_RTYP_TKR_TRAN 1
#define FES_RTYP_CAL_TRAN 2
#define FES_RTYP_ACD_TRAN 3
#define FES_DETC_TKR 0x00
#define FES_DETC_CAL 0x10
#define FES_DETC_ACD 0x20
#define FES_DTYP_SIMUL 1
#define FES_CAL_MRLENG 1400
#define FES_CAL_NRANGE 4
#define FES_TKR_MRLENG 14400 
#define FES_ACD_MRLENG 1000
#define FES_TKR_TRAN_SIZE 6


/* blw: to be deleted
class EbfContributorId
{
  public:
 
   enum 
   {
     TWR =  0,
     GEM = 16,
     ACD = 17
   };
};
*/




class FESOutput
{
  public:

//   inline FESOutput () { m_count = 0; return; }
FESOutput::FESOutput() :
    m_count(0),
    m_print(0),
    debug(false),
    firstEvt(true),
    m_GammaEvt(false)
{
    int i;
    for (i=0; i<16; ++i) {
        m_fpCAL[i] = 0;
        m_evtBufferCAL[i]=0;
        bLengthCAL[i]=0;
        }
    for (i=0; i<16; ++i)
        m_fpTKR[i] = 0;
        m_evtBufferTKR[i]=0;
        bLengthTKR[i]=0;
    for (i=0; i<4; ++i) {
        m_fpACD[i] = 0;
        m_evtBufferACD[i]=0;
        bLengthACD[i]=0;
        }
} 
         
   int          open   (const char    *fileName, const char *desc);

   unsigned int dumpTKR (const EbfTkrData *tkr, int nDeltaTime);
   unsigned int dumpCAL (const EbfCalData *cal, const Event::GltDigi &glt,int nDeltaTime);
   unsigned int dumpACD (const EbfAcdData *acd, int nDeltaTime);

   void         print  ()  const;
   int * writeTKR  (int tower, unsigned int *buff, int nBytes);
   int * writeCAL  (int tower, unsigned int *buff, int nBytes);
   int * writeACD  (int corner, unsigned int *buff, int nBytes);
   unsigned int  *writeFileHead(unsigned int *buff, int det, int tower, const char *desc);

   void  completeCAL(int nDelta);
   void  completeTKR(int nDelta);
   void  completeACD(int nDelta);
   void  countTrigger();
   void  dumpTriggerInfo();
   
   int fesFormatTKR(int timer, int nw_gtrc[8][9],
                                 int gtrc_addr[8][9],
                                 int tot_gtrc[8][9],
                                 int hit_gtrc[8][9][1536],
                                 int startLen,
                                 int twr);
   
   bool writeOutTKRTransitionVector(int detNumber, int num32, 
                                    int *word32,
                                    int *transVector);      

                                 
   unsigned int flush  ();
   unsigned int close  ();
   inline void setFirstEvtFlag(bool flg){ firstEvt=flg; return;}
   inline void setGammaFlag(bool flg)   { m_GammaEvt=flg; return;}
   inline void setVersion(int ver)
   {
      fesVersion= (ver==0 || ver==2) ? ver : 0; 
      return;
   }
   inline int   setPrint  (int printFlag) 
   {
       int old = m_print;
       m_print = printFlag;
       return old;
   }
             
   
   ~FESOutput ();
   
  /*
   | The unit of all size parameters is 32 - bit words
  */ 
  private:
    int                  m_count; /*!< Print flag                           */ 
    int                  m_print; /*!< Print flag                           */ 
    bool                   debug;    
    bool                 firstEvt;
    bool                 m_GammaEvt; /* Event is a Gamma Event */
    int                  fesVersion; /* Version 0: Cable ready, 2 Cable Ready w/ evtTim */

    int                  EvtTimebyteOffCAL; /* ptr to timing word for CAL data */
    int                  EvtTimebitOffCAL;  /* ptr to timing word for CAL data */
    int                  EvtTimebyteOffTKR; /* ptr to timing word for TKR data */
    int                  EvtTimebitOffTKR;  /* ptr to timing word for TKR data */
    int                  EvtTimebyteOffACD; /* ptr to timing word for ACD data */
    int                  EvtTimebitOffACD;  /* ptr to timing word for ACD data */

    int                  m_eventTrigger; /*Trigger Mask for event as determined by the FES */
    int                  m_eventCount[32]; /*Trigger Mask Counter*/

/* Files and counters for CAL file*/
    FILE               *m_fpCAL[16]; /*!< The file handle for CAL files    */   
    unsigned int    *m_evtBufferCAL[16]; /*!< The event buffer                     */
    unsigned int     bLengthCAL[16];
    unsigned int     m_CAL_evtCount[16];

/* Files and counters for TKR file*/
    FILE               *m_fpTKR[16]; /*!< The file handle for TKR files    */
    unsigned int    *m_evtBufferTKR[16]; /*!< The event buffer                     */
    unsigned int     bLengthTKR[16];
    unsigned int     m_TKR_evtCount[16];

/* Files and counters for ACD file*/
    FILE                *m_fpACD[4]; /*!< The file handle for ACD files    */
    unsigned int    *m_evtBufferACD[4]; /*!< The event buffer                     */
    unsigned int     bLengthACD[4];
    unsigned int     m_ACD_evtCount[4];

};




#endif




