#include <stdio.h>
#include <stdlib.h>

#include "EbfAcdData.h"
#include "EbfCalData.h"
#include "EbfTkrData.h"
#include "EbfGemData.h"
#include "FESOutput.h"

#undef DEBUG
#ifdef DEBUG
#define _DBG(statement) statement
#else
#define _DBG(statement)
#endif



    


static inline int addWord(int word, int wSize, unsigned int *buff, int lBuff,
                          int *byteOff, int *bitOff);

static inline int getWord(int *word, int wSize, unsigned int *buff, int lBuff,
                          int *byteOff, int *bitOff);
                   
static void printData (const unsigned int *beg, 
                       const unsigned int *end,
                       int rightMargin);



/**
 *
 *  @fn      FESOutput::open (const char    *fileName,
 *                            unsigned int maxEvtSize)
 *  @brief Initializes the FESOutput be opening the named file and allocating
 *         enough memory to hold at least 1 maximally sized event.
 *
 *  @param fileName   The name of the file to open.
 *  @param maxEvtSize The maximum size of a given event.
 *
 */
int FESOutput::open (const char     *fileName, const char *desc)
{

   debug = false;
   if (debug) printf("FESOutput::openning TKR file %s\n",fileName);


     
   /* Open the TKR output files */
   for(int tower=0; tower<16; tower++){

   /* initialize counter */
    m_TKR_evtCount[tower]=0;

   /* Initialize a chunk of memory to hold data for 1 tower */   
     m_evtBufferTKR[tower]  = (unsigned int *)malloc (2*FES_TKR_MRLENG);
     if (m_evtBufferTKR[tower] == 0) { return -1;  }


     char fileTKR[120];
     sprintf(fileTKR,"%s.tkr%i",fileName,tower);
     m_fpTKR[tower] = fopen (fileTKR, "wb");
     if (!m_fpTKR[tower])
     {
       /* Error in opening the output file.. */
       printf("FESOutput::ERROR Openning FES File for tower %i\n",tower);
       free (m_evtBufferTKR);
       return -2;
     }

/* Write TKR File Header */
     unsigned int *beg = m_evtBufferTKR[tower]; 
     m_evtBufferTKR[tower] = writeFileHead(m_evtBufferTKR[tower],1,tower,desc);
//     printData(beg,m_evtBufferTKR[tower],80);
     writeTKR(tower,beg,FES_FILEHEAD_SIZE);
     m_evtBufferTKR[tower] = beg;        
     
   }

   
   
   if (debug) printf("FESOutput::openning CAL file %s\n",fileName);

     
   /* Open the CAL output files */
   for(int tower=0; tower<16; tower++){


   /* initialize counter */
    m_CAL_evtCount[tower]=0;


   /* Initialize a chunk of memory to hold data for 1 tower */   
     m_evtBufferCAL[tower]  = (unsigned int *)malloc (2*FES_CAL_MRLENG);
     if (m_evtBufferCAL[tower] == 0) { return -1;  }

     char fileCAL[120];
     sprintf(fileCAL,"%s.cal%i",fileName,tower);
     m_fpCAL[tower] = fopen (fileCAL, "wb+");
     if (!m_fpCAL[tower])
     {
       /* Error in opening the output file.. */
       printf("ERROR Openning File for tower %i\n",tower);
       free (m_evtBufferCAL[tower]);
       return -2;
     }

/* Write CAL File Header */
     unsigned int *beg = m_evtBufferCAL[tower]; 
     m_evtBufferCAL[tower] = writeFileHead(m_evtBufferCAL[tower],2,tower,desc);
//     printData(beg,m_evtBufferCAL[tower],80);
     writeCAL(tower,beg,FES_FILEHEAD_SIZE);
     m_evtBufferCAL[tower] = beg;        
     
   }
   
   
   if (debug) printf("FESOutput::openning ACD file %s\n",fileName);

     
   /* Open the ACD output (4) files */
   for(int corner=0; corner<4; corner++){ 

   /* initialize counter */
     m_ACD_evtCount[corner]=0;

   /* Initialize a chunk of memory to hold data for 1 tower */   
     m_evtBufferACD[corner]  = (unsigned int *)malloc (2*FES_CAL_MRLENG);
     if (m_evtBufferACD[corner] == 0) { return -1;  }

     char fileACD[120];
     sprintf(fileACD,"%s.acd%i",fileName,corner);
     m_fpACD[corner] = fopen (fileACD, "wb");
     if (!m_fpACD[corner])
     {
       /* Error in opening the output file.. */
       printf("ERROR Openning File\n");
       free (m_evtBufferACD[corner]);
       return -2;
     }

/* Write ACD File Header */
     unsigned int *beg = m_evtBufferACD[corner]; 
     m_evtBufferACD[corner] = writeFileHead(m_evtBufferACD[corner],3,corner,desc);
//     printData(beg,m_evtBufferACD[corner],80);
     writeACD(corner,beg,FES_FILEHEAD_SIZE);
     m_evtBufferACD[corner] = beg;        
   }
     
   return 0;
}

unsigned int * FESOutput::writeFileHead(unsigned int *buff, int det, int tower, const char *desc)
/*
   DESCRIPTION
   -----------
   Composes the file header 

   PARAMETERS
   ----------


   RETURNS
   -------
   The compose header word.
*/
{

   unsigned int evtCount = 0;
   bool debug = false;
/* Form standard record header for File Head */
   unsigned int recHead = (FES_FILEHEAD_SIZE<<16) | (fesVersion<<8) | FES_RTYP_HEADER;
   if(debug) printf("Wd 0: %8.8x  %8.8x \n",(int)buff,recHead);
   *buff++ = recHead;

/* Form the detector Mask  (TKR or CAL)*/
   unsigned int detMask = 0;
   if(det!=3) detMask = (1 << tower);
   if(det==2) detMask = detMask << 16;  

/* Form the detector Mask for ACD */   
/*    --needs numbering convention **/
   unsigned int ACDdetMask = 0;
   if(det==3) ACDdetMask =1 << tower;

/* store lower 3 bytes of detector mask along with Data Type*/
   unsigned int wd = ((detMask & 0xffffff) << 8) | FES_DTYP_SIMUL;   
   if (debug) printf("Wd 1: %8.8x %8.8x \n",(int)buff,wd);
   *buff++ = wd;
   
/* store upper 1 byte of detector mask along with ACD Mask 
   and first two bytes of trans. count. */
   wd = ((evtCount*2 & 0xffff) << 16 ) | (ACDdetMask << 8) | ((detMask & 0xff000000)>>24);
   if (debug) printf("Wd 2: %8.8x %8.8x \n",(int)buff,wd);
   *buff++ = wd;
   
/* Store last two bytes of trans count and first two bytes of num events */
   wd = ((evtCount*2 & 0xffff0000) >> 16) | ((evtCount & 0xffff) << 16);
   if (debug) printf("Wd 3: %8.8x %8.8x \n",(int)buff,wd);
   *buff++ = wd;
   
/* Store last two bytes of event count and first two bytes of description field (blank)*/
   char descTot[80];
   sprintf(descTot,"%s",desc);
   for(int i=strlen(descTot); i<80; i++ ) strcat(descTot," ");
   wd = (descTot[1]<<24) | (descTot[0] << 16) | ((evtCount & 0xffff0000) >> 16);
   if (debug) printf("Wd 4: %8.8x %8.8x \n",(int)buff,wd);
   *buff++ = wd;
   
/* Store rest of description field plus padding (8 bytes) */
   for(int i=0; i<19; i++) *buff++ = descTot[4*i+2] | descTot[(4*i)+3]<<8 | descTot[(4*i)+4]<<16 | descTot[(4*i)+5]<<24 ;  
   *buff++ = 0xffff4321;   
      
   return (unsigned int *)buff;
}


/**
 *
 *  @fn   unsigned int FESOutput::dumpTKR (const EbfTkrData *tkr)
 *  @brief Converts the detector data blocks into EBF format.
 *  @param tkr  The TKR data block, contains all 16 towers.
 *  @return     ???
 *
 */
unsigned int FESOutput::dumpTKR(const EbfTkrData *tkr, int nDeltaTime){

   debug = false;
   if (debug) printf("FES dumpTKR Routine\n");

// Before dealing with this event complete the time stamp on the
// previous event and write it to the file
   if(!firstEvt) completeTKR(nDeltaTime);

   /*
    |
    | The numbering scheme is such that
    |
    |   X_0 lo =  0       Y_0 lo = 36
    |   X_0 hi =  1       Y_0 hi = 37
    |
    |   X_1 lo =  2       Y_1 lo = 38
    |   X_1 hi =  3       Y_1 hi = 39
    |
    |   X_2 lo =  4       Y_2 lo = 40
    |   X_2 hi =  5       Y_2 hi = 41
    |
    |   X_3 lo =  6       Y_3 lo = 42
    |   X_3 hi =  7       Y_3 hi = 43
    |
    |   X_4 lo =  8       Y_4 lo = 44
    |   X_4 hi =  9       Y_4 hi = 45
    |
    |   X_5 lo = 10       Y_5 lo = 46
    |   X_5 hi = 11       Y_5 hi = 47
    |
    |   X_6 lo = 12       Y_6 lo = 48
    |   X_6 hi = 13       Y_6 hi = 49
    |
    |   X_7 lo = 14       Y_7 lo = 50
    |   X_7 hi = 15       Y_7 hi = 51
    |
    |   X_8 lo = 16       Y_8 lo = 52
    |   X_8 hi = 17       Y_8 hi = 53
    |
    |   X_9 lo = 18       Y_9 lo = 54
    |   X_9 hi = 19       Y_9 hi = 55
    |
    |   X10 lo = 20       Y10 lo = 56
    |   X10 hi = 21       Y10 hi = 57
    |
    |   X11 lo = 22       Y11 lo = 58
    |   X11 hi = 23       Y11 hi = 59
    |
    |   X12 lo = 24       Y12 lo = 60
    |   X12 hi = 25       Y12 hi = 61
    |
    |   X13 lo = 26       Y13 lo = 62
    |   X13 hi = 27       Y13 hi = 63
    |
    |   X14 lo = 28       Y14 lo = 64
    |   X14 hi = 29       Y14 hi = 65
    |
    |   X15 lo = 30       Y15 lo = 66
    |   X15 hi = 31       Y15 hi = 67
    |
    |   X16 lo = 32       Y16 lo = 68
    |   X16 hi = 33       Y16 hi = 69
    |
    |   X17 lo = 34       Y17 lo = 70
    |   X17 hi = 35       Y17 hi = 71
   */

   /* 
    | Using the above information, fill a data structure giving the number
    | of the first layer on a controller. 
    /
    /  Note: The FES hardware assumes the cable order is 0 - 7. This is
    /        different than the ebf format output (and output of hardware),
    /        so we use a different "Begin" array.  
   */
   static const int Begin[EbfTkrTowerData::NumCables] = { 36,  37,  1,  0, 39,38, 2, 3};
                                                        

//   static const int Begin[EbfTkrTowerData::NumCables] = { 2,  0,  3,  1,
//                                                         38, 36, 39, 37  };
   
   unsigned int  accepts[EbfTkrTowerData::NumCables];
   int                                    nlayerEnds;
   int                                          boff;
   int           tots[EbfTkrTowerData::NumLayerEnds];


    /* Loop over the Tower Contributions */
   for (unsigned int itower = 0; itower < EbfCalData::NumTowers; itower++)
   {
      if (debug) printf("tower %d\n",itower);
      int nw_gtrcEmpty[8][9];
      int nw_gtrc[8][9];
      int tot_gtrc[8][9];
      int gtrc_addr[8][9];
      int hit_gtrc[8][9][1536];

      for (int ic=0; ic<8; ic++) {
         for (int ilay=0; ilay<9; ilay++) {
            nw_gtrc[ic][ilay] = 0;
         }
      }

      const EbfTkrTowerData *tkrTower = tkr->tower(itower);   

      /* 
       | First compute the accept bit list for the cable ends 
       | The output cable order is
       |      6. X odd  lo ( 2,  6, 10, 14, 18, 22, 26, 30, 34)
       |      3. X even lo ( 0,  4,  8, 12, 16, 20, 24, 28, 32)
       |      7. X odd  hi ( 3,  7, 11, 15, 19, 23, 27, 31, 35)
       |      2. X even hi ( 1,  5,  9, 13, 17, 21, 25, 29, 33)
       |      5. Y odd  lo (38, 42, 46, 50, 54, 58, 62, 66, 70)
       |      0. Y even lo (36, 40, 44, 48, 52, 56, 60, 64, 68)
       |      4. Y odd  hi (39, 43, 47, 51, 55, 59, 63, 67, 71)
       |      1. Y even hi (37, 41, 45, 49, 53, 57, 61, 65, 69)
      */
      nlayerEnds = 0;
      int cablePair = 0;
      for (unsigned int ocable = 0; ocable<EbfTkrTowerData::NumCables; ocable++)
      {
         if (debug) printf("cable %d\n",ocable);
          /*
           | mlayer is the accept mask for the current layer. The lowest layer
           | is placed in the MSB of the layer accept mask.
          */
          int mlayer       = (1 << (EbfTkrTowerData::NumLayerEndsPerCable - 1));
          int ilayerEndMin = Begin[ocable];
          int ilayerEndMax = ilayerEndMin 
                           + 4*EbfTkrTowerData::NumLayerEndsPerCable;
          int         list = 0;

          /* Pack the 9 layer accept bits for this cable */
          int layNum = 0;
	  
	  bool encode = false;
	  
          for (int ilayerEnd = ilayerEndMin;
               ilayerEnd     < ilayerEndMax;
               ilayerEnd    += 4, mlayer >>= 1)
          {
              nw_gtrc[ocable][layNum] = 0;
              nw_gtrcEmpty[ocable][layNum] = 0;
              gtrc_addr[ocable][layNum] = layNum;
              /* 
               | The number of layer ends with hits must be kept track of
               | so that one can pack one TOT per hit layer end. The tots
               | are sorted into an array so they can by jammed out in
               | the correct order.
              */
	      if(debug) printf("Checking LayerEnd Number %i which is layer %i on cable %i \n",ilayerEnd,layNum,ocable);
              if (tkrTower->m_nhits[ilayerEnd] || encode) 
              {
                  if (debug) printf("Hit found in tower %d cable %d [%d] tot %d hits %d\n",
                     itower,ocable,ilayerEnd,tkrTower->m_tots[ilayerEnd],
                     tkrTower->m_nhits[ilayerEnd]);
	       /* As far as I can tell these next three lines are not used BLW 17-Aug */	     
                  tots[nlayerEnds] = tkrTower->m_tots[ilayerEnd];
                  nlayerEnds      += 1;
                  list            |= mlayer;


                 /*
                  | Have some hits.
                  |  1. Compute a pointer to the strip ids on this layer.
                  |  2. Change the layer end into a just a layer number
                  |  3. Print the preamble for the hits
                 */
                 const unsigned short int *stripIds;
                 stripIds = tkrTower->m_data[ilayerEnd];       

                 /* Loop over all the strip hits on this layer end */
                 nw_gtrc[ocable][layNum] = tkrTower->m_nhits[ilayerEnd];
	         if(encode) nw_gtrc[ocable][layNum] = 64;
                 gtrc_addr[ocable][layNum] = layNum;
                 tot_gtrc[ocable][layNum] = tkrTower->m_tots[ilayerEnd];
		 /* if(encode) tot_gtrc[ocable][layNum] = */
                 for (int ihits = 0; ihits < nw_gtrc[ocable][layNum]; ihits++)
                 {
                     unsigned short int stripId = stripIds[ihits];
                     hit_gtrc[ocable][layNum][ihits] = stripId;
		     if(encode) hit_gtrc[ocable][layNum][ihits] = ( (ocable & 0x7) << 8 |
		     						    (layNum & 0xf) <<4  |
								    (ihits  & 0xf) ) ;
                     if (debug) printf("Tower %i Cable %i Layer %i hit %i   strip id %d\n",itower,ocable,layNum,ihits,hit_gtrc[ocable][layNum][ihits]);
                 }

              }
              layNum++;
          }

          accepts[ocable] = list;
          if (ocable%2 == 1) cablePair++;
          if (debug) printf("ocable %d cablepair %d\n",ocable,cablePair);
      }
   //
   // Now we have the data for this tower, let's format it
/*      int word32[1000];
      int word32Empty[1000];
      int transVector[6];
      int transVectorEmpty[6];
*/
      
// initialize the memory space
      size_t buffSize = 2*FES_TKR_MRLENG;  
      memset(m_evtBufferTKR[itower],0,buffSize);
      
      int timer = 0;       // appropriate time for real event
      int nBytes = fesFormatTKR(timer,nw_gtrc,gtrc_addr,tot_gtrc,hit_gtrc,0,itower);
      timer = 2*20;  // 2 x [20  x (50nsec)] = 2 microsec
      bLengthTKR[itower] = fesFormatTKR(timer,nw_gtrcEmpty,gtrc_addr,tot_gtrc,hit_gtrc,nBytes,itower);
 
//
// Write out the event Transition vector
/*      int detNumber = itower;
      if (writeOutTKRTransitionVector(detNumber,numWord32,word32,transVector)) {
         if(debug) printf("write was successful!\n");
      } else {
         printf("write failed!\n");
      }
//
// Write out the "empty" transition vector
      if (writeOutTKRTransitionVector(detNumber,numWord32Empty,word32Empty,transVectorEmpty)) {
         if(debug) printf("write was successful!\n");
      } else {
         printf("write failed!\n");
      }
*/      
   } 


return 0;
}


int FESOutput::fesFormatTKR(int timer,
                                 int nw_gtrc[8][9],
                                 int gtrc_addr[8][9],
                                 int tot_gtrc[8][9],
                                 int hit_gtrc[8][9][1536],
                                 int startLen,
                                 int twr
                                 )
{
   debug = false;
   int dsptr[8];
   int treq[8];
   int event[8];
   int nwordtot_11bit = 0;
	int nwordtot_15bit = 0;
	int nword_32bit = 0;

   int byteOff =startLen/4;
   int bitOff  =8*(startLen%4);
   
// Add blank words to hold the record header
   int HeadbyteOff =byteOff;
   int HeadbitOff  =bitOff;
   addWord(0,32,m_evtBufferTKR[twr],FES_TKR_MRLENG,&byteOff,&bitOff);
   addWord(0,8,m_evtBufferTKR[twr],FES_TKR_MRLENG,&byteOff,&bitOff);

// If this is version 2 (evt timing) add a blank 32 bit word for the timing
// if second transition record for event, put in timer.
   if(fesVersion==2) {
     if(startLen==0) {
       EvtTimebyteOffTKR = byteOff;
       EvtTimebitOffTKR  = bitOff;
       addWord(0,32,m_evtBufferTKR[twr],FES_TKR_MRLENG,&byteOff,&bitOff);
     } else {
       addWord(timer,32,m_evtBufferTKR[twr],FES_TKR_MRLENG,&byteOff,&bitOff);
     }  
   }

// Add blank words for the transition vectors
   int TranVbyteOff = byteOff;
   int TranVbitOff  = bitOff;
   for(int i=0; i<FES_TKR_TRAN_SIZE; i++) addWord(0,32,m_evtBufferTKR[twr],FES_TKR_MRLENG,&byteOff,&bitOff);  


   for (int ic=0; ic<8; ic++) {
      if (debug) printf("cable %d\n",ic);
      event[ic] = 0;
      treq[ic] = 0;
//
// Make a list of 11 bit words for this cable
      int nword_11bit = 0;
      int word_11bit[1000];
      for (int ilay=0; ilay<9; ilay++) {
         if (debug) printf("Layer %d hits %d\n",ilay,nw_gtrc[ic][ilay]);
         int header = gtrc_addr[ic][ilay];
         if (nw_gtrc[ic][ilay] > 0) header += 16;
         if (debug) printf("header 0x%x gtrc_addr %d\n",header,gtrc_addr[ic][ilay]);
         word_11bit[nword_11bit] = header;
         nword_11bit++;
         if (nw_gtrc[ic][ilay] > 0) {
            event[ic] = 1;
            treq[ic] = treq[ic] | (1<<ilay);
            word_11bit[nword_11bit] = nw_gtrc[ic][ilay] + 1;
            nword_11bit++;
            for (int ihits = 0; ihits < nw_gtrc[ic][ilay]; ihits++)
            {
               if (debug) printf("   hit %d %d\n",ihits,hit_gtrc[ic][ilay][ihits]);
               word_11bit[nword_11bit] = hit_gtrc[ic][ilay][ihits];
               nword_11bit++;
            }
            if (debug) printf("    TOT %d %x\n",tot_gtrc[ic][ilay],tot_gtrc[ic][ilay]);
            word_11bit[nword_11bit] = tot_gtrc[ic][ilay];
            nword_11bit++;
            
         }
      }
      nwordtot_11bit += nword_11bit;
      if (debug) printf("\nTotal number of 11 bit words %d\n",nword_11bit);
      for (int i11=0; i11<nword_11bit; i11++) {
         if (debug) printf("   i11 %d iword %x\n",i11,word_11bit[i11]);
      }
//
// Loop over the 11 bit words, and merge into 15 bit words
      int ibitStart = 0;
      int ibitEnd = 14;
      int ibitCurrent = 0;
      int nword_15bit = 0;
      int word15[1000];
      word15[0] = 0;
      bool lastWord = false;
      int current_nword_15bit = nword_15bit;
      nword_15bit++;
      word15[nword_15bit]=0;
      int ibit15=0;
      for (int i11=0; i11<nword_11bit; i11++) {
         if (debug) printf("11bit word %x\n",word_11bit[i11]);
         for (int ibit=0; ibit<11; ibit++) {
            if (word_11bit[i11] & (1<< ibit)) {
               word15[nword_15bit] = word15[nword_15bit] | (1<<ibit15);
            }
            if (debug) printf("bit count 11b %d 15b %d word15 %x\n",ibit,ibit15,word15[nword_15bit]);
            ibit15++;
            if (ibit15 > 14) {
               ibit15=0;
               nword_15bit++;
               word15[nword_15bit]=0;
            }
         }
      }
      nword_15bit++;
      if (debug) {
        printf("Total number of 15 bit words %d\n",nword_15bit);
        for (int i15=0; i15<nword_15bit; i15++) {
          printf("   i15 %d iword %x\n",i15,word15[i15]);
        }
      }  
//
// If there is an odd number of words, then add one and pad with zeroes
      if (nword_15bit%8 != 0) {
         int jpad = 8 - nword_15bit%8;
         if (debug) printf("Number of 15 bit words: %d, pad to %d\n",nword_15bit,jpad);
	      if (jpad != 0) {
            for (int ipad=0; ipad<jpad; ipad++) {
		         word15[nword_15bit] = 0;
		         nword_15bit++;
            }
	      }
      }
      int num8words = nword_15bit/8;
      word15[current_nword_15bit] = num8words-1;
      if (debug) printf("contents of length word array %d; address %d\n",
               word15[current_nword_15bit],current_nword_15bit);
      nwordtot_15bit += nword_15bit;
      if (debug) printf("Total number of 15 bit words %d\n",nword_15bit);
      for (int i15=0; i15<nword_15bit; i15++) {
         if (debug) printf("   i15 %d iword %x\n",i15,word15[i15]);
      }
//
// The current value of the 32 bit word counter is the pointer to the event data
      dsptr[ic] = num8words;
//
// Now merge the 15 bit words into 32 bit words (bits 15 and 31 set to zero)
//      int current_nword_32bit = nword_32bit;
//      nword_32bit++;
	   for (int iword15=0; iword15<nword_15bit; iword15=iword15+2) {
		   int word = 0;
         if (debug) printf("iword %d word15[0] %x word15[1] %x\n",
            iword15,word15[iword15],word15[iword15+1]);
		   word = (word15[iword15+1] << 16) | word15[iword15];
         if (debug) printf("word %x\n",word);
// Leave Space for the header and the transition vectors         
//		   word32[nword_32bit] = word;
         addWord(word,32,m_evtBufferTKR[twr],FES_TKR_MRLENG,&byteOff,&bitOff);         
		   nword_32bit++;
	   }
//      word32[current_nword_32bit] = nword_32bit - current_nword_32bit;
   }

   if (debug) printf("num11bit %d num15bit %d num32bit %d\n",
      nwordtot_11bit,nwordtot_15bit,nword_32bit);
 
// Get the length of this record 
   int recLength = byteOff*4 + bitOff/8;

//
// Now for the transition vectors
   int eventPair0=0;
   int eventPair1=0;
   int eventPair2=0;
   int eventPair3=0;
//
// 1st/2nd cable pair
   if ((event[0] == 1) || (event[1] == 1) ||
      (event[2] == 1) || (event[3] == 1) ||
      (event[4] == 1) || (event[5] == 1) ||
      (event[6] == 1) || (event[7] == 1)) {
      eventPair0 = 1;
      eventPair1 = 1;
      eventPair2 = 1;
      eventPair3 = 1;
   }
   
   int transVector[6];
   transVector[0] = (dsptr[1]<<25) | (treq[1]<<16) | (dsptr[0]<<9) | treq[0];
   transVector[1] = (dsptr[2]<<25) | (treq[2]<<16) | (eventPair0<<14) | timer;
   transVector[2] = (eventPair1<<30) | (timer<<16) | (dsptr[3]<<9) | treq[3];
   if (debug) printf("TREQ0: %x  DATP0: %x  TREQ1 %x  DATP1: %x\n",
         treq[0],dsptr[0],treq[1],dsptr[1]);
   if (debug) printf("Time0: %x  Event0: %x\n",eventPair0,timer);
   if (debug) printf("TREQ2: %x  DATP2: %x  TREQ3 %x  DATP3: %x\n",
         treq[2],dsptr[2],treq[3],dsptr[3]);
   if (debug) printf("Time1: %x  Event1: %x\n",eventPair1,timer);
//
// 3rd/4th cable pair
   transVector[3] = (dsptr[5]<<25) | (treq[5]<<16) | (dsptr[4]<<9) | treq[4];
   transVector[4] = (dsptr[6]<<25) | (treq[6]<<16) | (eventPair2<<14) | timer;
   transVector[5] = (eventPair3<<30) | (timer<<16) | (dsptr[7]<<9) | treq[7];
   if (debug) printf("TREQ4: %x  DATP4: %x  TREQ5 %x  DATP5: %x\n",
         treq[4],dsptr[4],treq[5],dsptr[5]);
   if (debug) printf("Time2: %x  Event2: %x\n",eventPair2,timer);
   if (debug) printf("TREQ6: %x  DATP6: %x  TREQ7 %x  DATP7: %x\n",
         treq[6],dsptr[6],treq[7],dsptr[7]);
   if (debug) printf("Time3: %x  Event3: %x\n",eventPair3,timer);

// Attach the record header
   addWord(FES_RTYP_TKR_TRAN,8,m_evtBufferTKR[twr],FES_TKR_MRLENG,&HeadbyteOff,&HeadbitOff);
   addWord(fesVersion,8,m_evtBufferTKR[twr],FES_TKR_MRLENG,&HeadbyteOff,&HeadbitOff);
   addWord(recLength-startLen,16,m_evtBufferTKR[twr],FES_TKR_MRLENG,&HeadbyteOff,&HeadbitOff);

// add the detector number
   int detID = twr | FES_DETC_TKR;
   if(m_GammaEvt) detID |= 0x80;
   addWord(detID,8,m_evtBufferTKR[twr],FES_TKR_MRLENG,&HeadbyteOff,&HeadbitOff);   


// Put the Transition Vector at the beginning of the record
   for(int i=0; i<FES_TKR_TRAN_SIZE; i++) addWord(transVector[i],32,m_evtBufferTKR[twr],FES_TKR_MRLENG,&TranVbyteOff,&TranVbitOff);   

   return recLength;
}


bool FESOutput::writeOutTKRTransitionVector(int detNumber,
                                          int num32, 
                                          int *word32,
                                          int *transVector)      
{

   int recordType = 1;
   int formatVersion = 0;
   int recordLength = 5 + 4*6+4*num32;
   fwrite(&recordType, 1, 1, m_fpTKR[detNumber]);
   fwrite(&formatVersion, 1 , 1, m_fpTKR[detNumber]);
   fwrite(&recordLength, 1, 2, m_fpTKR[detNumber]);
   fwrite(&detNumber, 1, 1, m_fpTKR[detNumber]);
   fwrite(transVector, 1, 4*6, m_fpTKR[detNumber]);
   fwrite(word32, 1, 4*num32, m_fpTKR[detNumber]);

   
   return true;
}

/**
 *
 *  @fn   unsigned int FESOutput::completeTKR (int nDelta)
 *  @brief finishes write out of CAL FES data
 *  @param nDelta  number of clock ticks since last event
 *  @return     void
 *
 */
void FESOutput::completeTKR(int nDeltaTime){

// Loop over towers
  for(int tower=0; tower<16; tower++){


// If this is version 2, push in the delta time (LAT clock ticks to next event)
     if(fesVersion==2) {
         unsigned int timeWd = (1 << 31) | (nDeltaTime & 0x3fffffff);
//         printf("Event Time Word 0x%8.8x , EvtTimebyteOffCAL %i  EvtTimebitOff %i \n",timeWd,EvtTimebyteOffCAL,EvtTimebitOffCAL);
         int byteOff = EvtTimebyteOffTKR;
         int bitOff = EvtTimebitOffTKR;
         addWord(timeWd,32,m_evtBufferTKR[tower],FES_TKR_MRLENG,&byteOff,&bitOff);
     } 
  
     writeTKR(tower,m_evtBufferTKR[tower],bLengthTKR[tower]);

// Keep the number of events written into this file. At the end of the job
// we will rewind the file and put the info in the header.
     m_TKR_evtCount[tower]++;     

     
  }

  return;
}

/**
 *
 *  @fn   unsigned int FESOutput::dumpCAL (const EbfCalData *cal)
 *  @brief Converts the detector data blocks into EBF format.
 *  @param cal  The CAL data block, contains all 16 towers.
 *  @return     ???
 *
 */
unsigned int FESOutput::dumpCAL(const EbfCalData *cal, const Event::GltDigi &glt, int nDeltaTime){


// Before dealing with this event complete the time stamp on the
// previous event and write it to the file
   if(!firstEvt) completeCAL(nDeltaTime);

   bool debug = false;
   bool encode = false;
   if (debug) printf("FES dumpCAL Routine\n");

   for(int tower=0; tower<16; tower++) {

//      debug = tower==7 ? true : false;

      if(debug) printf("  **** Tower %i ******\n\n",tower);

/* Calorimeter Tower */
/* ------------------*/   
      const EbfCalTowerData *calTower = cal->tower(tower);

      size_t buffSize = 2*FES_CAL_MRLENG;  
      memset(m_evtBufferCAL[tower],0,buffSize);

/* Offsets used for adding words */         
      int byteOff = 0;
      int bitOff = 0;

/* Start buffer for transition record...leave room at top for record header 5 bytes word */
      unsigned int *beg = m_evtBufferCAL[tower];
      int zero = 0;
      int HeaderByteOff = byteOff;
      int HeaderBitOff = bitOff;
      addWord(zero,32,m_evtBufferCAL[tower],FES_CAL_MRLENG,&byteOff, &bitOff);
      addWord(zero,8,m_evtBufferCAL[tower],FES_CAL_MRLENG,&byteOff, &bitOff);

/* If event time version (2) we need a 4 byte word for the time stamp */
      EvtTimebyteOffCAL = byteOff;
      EvtTimebitOffCAL = bitOff;
      if(fesVersion==2) addWord(zero,32,m_evtBufferCAL[tower],FES_CAL_MRLENG,&byteOff, &bitOff);

/* Add the data Must leave room for the transition records 4-32bit words = 16 Bytes */
      int TVecByteOff = byteOff;
      int TVecBitOff  = bitOff;  
      addWord(zero,32,m_evtBufferCAL[tower],FES_CAL_MRLENG,&byteOff, &bitOff);
      addWord(zero,32,m_evtBufferCAL[tower],FES_CAL_MRLENG,&byteOff, &bitOff);
      addWord(zero,32,m_evtBufferCAL[tower],FES_CAL_MRLENG,&byteOff, &bitOff);
      addWord(zero,32,m_evtBufferCAL[tower],FES_CAL_MRLENG,&byteOff, &bitOff);   


/* Expand adc hits out into an array without sparsification
 * and arrange data a long the lines of FES Cables and Layers
 *          OK Ebf stores the adc values packed in a single 32 bit word
 *               bits 31-28 (4): log number 0-11
 *               bits 27-26 (2): Range for -end
 *               bits 25-14 (12): Adc Value for -end
 *               bits 13-12 (2): Range for +end
 *               bits 11-0  (12): Adc Value for +end
 * 
 *   adc[4][4][4][12]
 *     4 FES cables:
 *     4 layers (FES #): 
 *     4 ranges: Note Ebf only stores a single range...
 *     12 logs:
 *
 *   dav[4][4]:
 *     4 FES cables:
 *     4 layers:
 *
 *   acc[4][4]:
 *     4 FES cables:
 *     4 layers:
 *
 *   rngLsb[4][4]: 12 bits for least sign. bit of range
 *     4 FES cables:
 *     4 layers:
 *
 *   rngMsb[4][4]: 12 bits for most sign. bit of range
 *     4 FES cables:
 *     4 layers:
 *
 *   Ebf to FES Cable/Layer mapping:
 *      Efb Layer:       FES Cable/Layer
 *         0 (+x)            0      0
 *         0 (-x)            2      0
 *         2 (+x)            0      1
 *         2 (-x)            2      1
 *         4 (+x)            0      2
 *         4 (-x)            2      2
 *         6 (+x)            0      3
 *         6 (-x)            2      3
 *
 *         1 (+y)            1      0
 *         1 (-y)            3      0
 *         3 (+y)            1      1
 *         3 (-y)            3      1
 *         5 (+y)            1      2
 *         5 (-y)            3      2
 *         7 (+y)            1      3
 *         7 (-y)            3      3  

*/
      unsigned int adc[4][4][4][12];
      memset(adc, 0x0, sizeof adc);

      unsigned int dav[4][4];
      memset(dav, 0x0, sizeof dav);

      unsigned int acc[4][4];
      memset(acc, 0x0, sizeof acc);

      unsigned int rngLsb[4][4];
      memset(rngLsb, 0x0, sizeof rngLsb);

      unsigned int rngMsb[4][4];
      memset(rngMsb, 0x0, sizeof rngMsb);

      unsigned int treq[4];
      memset(treq, 0x0, sizeof treq);

      unsigned int nhit[4][4][2];
      memset(nhit, 0x0, sizeof nhit);
      
      unsigned int layerHits[4][4][6][2];
      memset(layerHits, 0x0, sizeof nhit);

// Initial the layer hit arrays and counters.
      for(int FESCable=0; FESCable<4; FESCable++) {
          for(int FESlayer=0; FESlayer<4; FESlayer++) {
             for(int hit=0; hit<6; hit++){
                layerHits[FESCable][FESlayer][hit][0]=0xF;
                layerHits[FESCable][FESlayer][hit][1]=0xF;                 
             }
             nhit[FESCable][FESlayer][0] = 0;
             nhit[FESCable][FESlayer][1] = 0;
          }
      }

      for(int layer=0; layer<calTower->NumLayers; layer++) {
         int numHits = calTower->m_layers[layer].m_nhits;
         unsigned int layerHitMask = calTower->m_layers[layer].m_msk;
         if(encode) layerHitMask = 0xfff;
         if(debug) printf("Number of Hits for Layer %i is %i Hit mask 0x%3.3x\n",layer,calTower->m_layers[layer].m_nhits,layerHitMask);
         for(int nlog=0; nlog<12; nlog++) {
// If no hit in this layer try the next log         
             if((layerHitMask & (1 << nlog)) == 0) continue;

// If look for all four ranges if necesary
             int lastRange = (cal->FourRangeReadOut()) ? 4:1;
             for(int readoutRange=0; readoutRange<lastRange; readoutRange++) {
               unsigned int adcWord = calTower->m_layers[layer].m_xtals[nlog][readoutRange];             
             
               int log = (adcWord & (0xf << 28))>>28;
               int rangeN = (adcWord & (0x3 << 26))>>26;
               int adcN = (adcWord & (0xfff << 14))>>14;          
               int rangeP = (adcWord & (0x3 << 12))>>12;
               int adcP = (adcWord & 0xfff);

//               if(debug) printf("Unpacked Words Tower %i ADC 0x%8.8x  : log %i  rangeP %i  adcP 0x%8.8x rangeN %i adcN 0x%8.8x \n",tower,adcWord,log,rangeP,adcP,rangeN,adcN);

/* FESCable is 0 (+x) or 1 (+y)...-x and -y add 2
  *  FESLayer is 0-3 divide Ebf layer by 2 (integer math) */
               int FESCable = layer%2;
               int FESlayer = layer/2;

               adc[FESCable][FESlayer][rangeP][log] = adcP;
               adc[FESCable+2][FESlayer][rangeN][log] = adcN;     

/* Certian things only once for the four ranges */
               if(readoutRange==0) {


/* Keep track of the number of hits within a half layer and which log they are on POSITIVE END*/
                 int halfLayer = log < 6 ? 0 : 1;
                 layerHits[FESCable][FESlayer][nhit[FESCable][FESlayer][halfLayer]][halfLayer] = halfLayer==1 ? log-6 : log;
                 nhit[FESCable][FESlayer][halfLayer]++;

/* Keep track of the number of hits within a half layer and which log they are on NEGATIVE END*/
                 halfLayer = log < 6 ? 0 : 1;
                 layerHits[FESCable+2][FESlayer][nhit[FESCable+2][FESlayer][halfLayer]][halfLayer] = halfLayer==1 ? log-6 : log;
                 nhit[FESCable+2][FESlayer][halfLayer]++;

/* Do we need to add something for the 4 range read out? */            

/* Keep track of dav, acc, and rng bits */
//                 debug = tower==7 ? true : false;
                 if(adcP != 0 || encode) {
                    dav[FESCable][FESlayer] |= (0x1 << log);
                    acc[FESCable][FESlayer] |= (0x1 << log);
                    rngLsb[FESCable][FESlayer] |= ( (rangeP & 0x1) << log);
                    rngMsb[FESCable][FESlayer] |= ( ((rangeP & 0x2)>>1) << log);
                    if (debug) printf("Hits Found: Layer %i Log %i FESCable %i FESLayer %i RangeP %i ADCP 0x%8.8x RangeN %i ADCN 0x%8.8x MSB 0x%8.8x LSB 0x%8.8x\n",layer,log,FESCable,FESlayer,rangeP,adcP,rangeN,adcN,rngMsb[FESCable][FESlayer],rngLsb[FESCable][FESlayer]);

// Check trigger thresholds
                    if(glt.getCALLOtrigger(idents::CalXtalId(tower,layer,log,0)))  treq[FESCable] |= (1 << FESlayer);
                    if(glt.getCALHItrigger(idents::CalXtalId(tower,layer,log,0)))  treq[FESCable] |= (1 << (FESlayer+4));
                 }
                 if(adcN != 0 || encode) {
                    dav[FESCable+2][FESlayer] |= (0x1 << log);
                    acc[FESCable+2][FESlayer] |= (0x1 << log);
                    rngLsb[FESCable+2][FESlayer] |= ( (rangeN & 0x1) << log);
                    rngMsb[FESCable+2][FESlayer] |= ( ((rangeN & 0x2)>>1) << log);
                    if (debug && adcP == 0) printf("Hits Found: Layer %i Log %i FESCable %i FESLayer %i RangeP %i ADCP 0x%8.8x RangeN %i ADCN 0x%8.8x MSB 0x%8.8x LSB 0x%8.8x \n",layer,log,FESCable,FESlayer,rangeP,adcP,rangeN,adcN,rngMsb[FESCable][FESlayer],rngLsb[FESCable][FESlayer]);

// Check trigger thresholds
                    if(glt.getCALLOtrigger(idents::CalXtalId(tower,layer,log,1)))  treq[FESCable+2] |= (1 << FESlayer);
                    if(glt.getCALHItrigger(idents::CalXtalId(tower,layer,log,1)))  treq[FESCable+2] |= (1 << (FESlayer+4));                
                 } 
//                 debug = false;         
               }// readoutRange==0
            }//loop over ranges   
         }
      }

/*
* Loop over cables and fill the data buffer
*/
      for(int FESCable=0; FESCable<4; FESCable++) {
         if (debug) printf("CABLE: %i\n",FESCable);
/* 
*  Header for cable with dav bits, acc bits, and range bits
*/
         int dav_0_2a = ((dav[FESCable][2] & 0xff)  << 24)| 
                         ((dav[FESCable][1] & 0xfff) << 12)| 
                          (dav[FESCable][0] & 0xfff);
         if (debug) printf("Sum Word 0: 0x%8.8x\n",dav_0_2a);
         addWord(dav_0_2a,32,m_evtBufferCAL[tower],FES_CAL_MRLENG,&byteOff, &bitOff);

         int dav_2b_3_acc_0_1a =((acc[FESCable][1] & 0xf)   << 28)| 
                                ((acc[FESCable][0] & 0xfff) << 16)| 
                                ((dav[FESCable][3] & 0xfff) << 4 )|
                                ((dav[FESCable][2] & 0xf00) >> 8 ); 
         if (debug) printf("Sum Word 1: 0x%8.8x\n",dav_2b_3_acc_0_1a);
         addWord(dav_2b_3_acc_0_1a,32,m_evtBufferCAL[tower],FES_CAL_MRLENG,&byteOff, &bitOff);

         int acc_1b_3 = ( acc[FESCable][3] << 20)| 
                        ((acc[FESCable][2] & 0xfff) << 8)| 
                        ((acc[FESCable][1] & 0xff0) >> 4 );   
         if (debug) printf("Sum Word 2: 0x%8.8x\n",acc_1b_3);
         addWord(acc_1b_3,32,m_evtBufferCAL[tower],FES_CAL_MRLENG,&byteOff, &bitOff);      

         int rLsb_0_2a = ((rngLsb[FESCable][2] & 0xff)  << 24)| 
                         ((rngLsb[FESCable][1] & 0xfff) << 12)| 
                          (rngLsb[FESCable][0] & 0xfff);
         if (debug) printf("Sum Word 3: 0x%8.8x\n",rLsb_0_2a);
         addWord(rLsb_0_2a,32,m_evtBufferCAL[tower],FES_CAL_MRLENG,&byteOff, &bitOff);

         int rLsb_2b_3_rMsb_0_1a = ((rngMsb[FESCable][1] & 0xf)   << 28)| 
                                   ((rngMsb[FESCable][0] & 0xfff) << 16)| 
                                   ((rngLsb[FESCable][3] & 0xfff) << 4 )|
                                   ((rngLsb[FESCable][2] & 0xf00) >> 8 ); 
         if (debug) printf("Sum Word 4: 0x%8.8x\n",rLsb_2b_3_rMsb_0_1a); 
         addWord(rLsb_2b_3_rMsb_0_1a,32,m_evtBufferCAL[tower],FES_CAL_MRLENG,&byteOff, &bitOff);

         int rMsb_1b_3 = ( rngMsb[FESCable][3] << 20)| 
                         ((rngMsb[FESCable][2] & 0xfff) << 8)| 
                         ((rngMsb[FESCable][1] & 0xff0) >> 4 );
         if (debug) printf("Sum Word 5: 0x%8.8x at nByte %i and nBit %i\n",rMsb_1b_3,byteOff,bitOff);
         addWord(rMsb_1b_3,32,m_evtBufferCAL[tower],FES_CAL_MRLENG,&byteOff, &bitOff);       

/* Get the data available summary */
         unsigned int davSum = 0;
         for(int FESlayer=0; FESlayer<4; FESlayer++){
             davSum |= ( (dav[FESCable][FESlayer] & 0x3F) | (dav[FESCable][FESlayer] & 0xFC0)>>6 ); 
         }
// Count the number of 48 bit words for FES memory...start with 4
// because of the cable summary words.      
/*         int n48bitWds = 4;
//      printf("davSum 0x%8.8x\n",davSum);
         for(int range=0; range<FES_CAL_NRANGE; range++){
             for(int shift=6; shift>=0; shift -= 6) {
                for(int halfLog=0; halfLog<6; halfLog++){
//                 printf("Bit Check: %i\n",(davSum & (0x1 << halfLog)) >> halfLog);
                    if( (davSum & (0x1 << halfLog)) >> halfLog ) {
                        for(int FESlayer=0; FESlayer<4; FESlayer++) {
                           int cRange = (((rngMsb[FESCable][FESlayer]>>halfLog)<<1 | (rngLsb[FESCable][FESlayer]>>halfLog))&0x3) + range;
                           if(cRange >= FES_CAL_NRANGE) cRange -= FES_CAL_NRANGE;
                           unsigned int adcWrdA = (adc[FESCable][FESlayer][cRange][halfLog]>>shift)&0x3f;
                           addWord(adcWrdA,6,m_evtBufferCAL[tower],FES_CAL_MRLENG,&byteOff, &bitOff);
                           unsigned int adcWrdB = (adc[FESCable][FESlayer][cRange][halfLog+6]>>shift)&0x3f;
                           addWord(adcWrdB,6,m_evtBufferCAL[tower],FES_CAL_MRLENG,&byteOff, &bitOff);
                           if (debug) printf("Layer %i HalfLog %i Shift %i Range %i  cRange %i Wrd A 0x%8.8x Wrd B 0x%8.8x\n",
                                   FESlayer,halfLog,shift,range,cRange,adcWrdA,adcWrdB);                                
                        } 
                        n48bitWds++;
                    } 
                }
             }
         }
*/

// Get the maximum number of hits on a layer for this cable
         int maxCableHit = 0;
         for(int FESlayer=0; FESlayer<4; FESlayer++) {
//             if(debug) for(int n=0; n<6; n++) printf("C %i L %i Hit %i Hit0 %i Hit1 %i \n",FESCable,FESlayer,nhit,
//                    layerHits[FESCable][FESlayer][n][0],layerHits[FESCable][FESlayer][n][1]);
             int nhit = 5;
//             if(debug) printf("Cable %i Layer %i hit %i Hits0 %i Hits1 %i\n",FESCable,FESlayer,nhit,
//                    layerHits[FESCable][FESlayer][nhit][0],layerHits[FESCable][FESlayer][nhit][1]);
             while ( (layerHits[FESCable][FESlayer][nhit][0] == 0xF) & 
                     (layerHits[FESCable][FESlayer][nhit][1] == 0xF) &
                     (nhit>-1) ) {
                    nhit--;
//             if(debug) printf("Cable %i Layer %i hit %i Hits0 %i Hits1 %i\n",FESCable,FESlayer,nhit,
//                    layerHits[FESCable][FESlayer][nhit][0],layerHits[FESCable][FESlayer][nhit][1]);
             }       
             if(nhit+1 > maxCableHit) maxCableHit = nhit+1;
         }
         if(debug) printf("Maximum Hits on Cable %i is %i\n",FESCable,maxCableHit);

// Count the number of 48 bit words for FES memory...start with 4
// because of the cable summary words.      
         int n48bitWds = 4;
         for(int range=0; range<FES_CAL_NRANGE; range++){
             for(int shift=6; shift>=0; shift -= 6) {
                for(int nhit=0; nhit<maxCableHit; nhit++){
                   for(int FESlayer=0; FESlayer<4; FESlayer++){

// Where is the first hit? (Note: if there is no hit layerHits[][][][] = 0xF)
//                      int halfLog1 = layerHits[FESCable][FESlayer][nhit][0];
//                      int halfLog2 = layerHits[FESCable][FESlayer][nhit][1];
//                      int halfLog  = halfLog1 < halfLog2 ? halfLog1 : halfLog2;

// First get the range of the hit (Is this correct?)
//                      int cRange = (((rngMsb[FESCable][FESlayer]>>halfLog)<<1 | (rngLsb[FESCable][FESlayer]>>halfLog))&0x3) + range;
//                      if(cRange >= FES_CAL_NRANGE) cRange -= FES_CAL_NRANGE;

// Where is the first hit? (Note: if there is no hit layerHits[][][][] = 0xF)
                      int halfLog1 = layerHits[FESCable][FESlayer][nhit][0];
                      int halfLog2 = layerHits[FESCable][FESlayer][nhit][1];

// First get the range of the hit 
                      int cRange1 = (( ((rngMsb[FESCable][FESlayer]>>halfLog1)&0x1)<<1 | ((rngLsb[FESCable][FESlayer]>>halfLog1)&0x1))&0x3) + range;
                      if(cRange1 >= FES_CAL_NRANGE) cRange1 -= FES_CAL_NRANGE;
//                      if(debug) printf("Msb %i LsB %i \n",(rngMsb[FESCable][FESlayer]>>halfLog1)&0x1,(rngLsb[FESCable][FESlayer]>>halfLog1)&0x1);

// First get the range of the hit 
                      int cRange2 = (( ((rngMsb[FESCable][FESlayer]>>(halfLog2+6))&0x1)<<1 | ((rngLsb[FESCable][FESlayer]>>(halfLog2+6))&0x1))&0x3) + range;
                      if(cRange2 >= FES_CAL_NRANGE) cRange2 -= FES_CAL_NRANGE;
//                      if(debug) printf("Msb %i LsB %i \n",(rngMsb[FESCable][FESlayer]>>(halfLog2+6))&0x1,(rngLsb[FESCable][FESlayer]>>(halfLog2+6))&0x1);


                      unsigned int adcWrdA = halfLog1 < 6 ? (adc[FESCable][FESlayer][cRange1][halfLog1]>>shift)&0x3f : 0;
                      addWord(adcWrdA,6,m_evtBufferCAL[tower],FES_CAL_MRLENG,&byteOff, &bitOff);
                      unsigned int adcWrdB = halfLog2 < 6 ? (adc[FESCable][FESlayer][cRange2][halfLog2+6]>>shift)&0x3f : 0;
                      addWord(adcWrdB,6,m_evtBufferCAL[tower],FES_CAL_MRLENG,&byteOff, &bitOff);
                      if (debug) printf("Layer %i HalfLog1 %i  HalfLog2 %i  Shift %i Range %i  cRange1 %i  Wrd A 0x%8.8x cRange2 %iWrd B 0x%8.8x\n",
                              FESlayer,halfLog1,halfLog2,shift,range,cRange1,adcWrdA,cRange2,adcWrdB);                                

                                              
                   }
                  n48bitWds++;  
                }
             }
         }


/* Pad out the data...to give an integer number of 8-48 bit words...we need 4x48 bit word = 8x24 bit words */
         for(int i=0; i<8; i++) addWord(0,24,m_evtBufferCAL[tower],FES_CAL_MRLENG,&byteOff, &bitOff);
/* Increment the number of 48 bit words for this padding */
         n48bitWds += 4;   

/* Add the transition data */
         int datpt  = n48bitWds/8;
         if (debug) printf("n48bitWds %i  datpt %i\n",n48bitWds,datpt);
         int time   = 0;      
         unsigned int tVec = ( (1<<30) | (time & 0x3fff)<<16) | ((datpt & 0x7)<<8) | (treq[FESCable] & 0xff);
         addWord(tVec,32,m_evtBufferCAL[tower],FES_CAL_MRLENG,&TVecByteOff, &TVecBitOff);
//      printf("tvec 0x%8.8x\n",tVec);  

      } //end loop over FES Cable 

/* Go back and add the record length now that we know it */   
//   unsigned int lengthBuff = (unsigned int *)m_evtBufferCAL[tower] - (unsigned int *)beg;
      unsigned int nBytes = byteOff*4 + bitOff/8;
      unsigned int recHead = ((nBytes)<<16) | (fesVersion<<8) | FES_RTYP_CAL_TRAN;
      if (debug) printf("record header 0x%8.8x byteOff %i bitOff %i nBytes %i\n",recHead,byteOff,bitOff,nBytes);

      addWord(recHead,32,m_evtBufferCAL[tower],FES_CAL_MRLENG,&HeaderByteOff, &HeaderBitOff);
      int detID =  tower | FES_DETC_CAL;
      if(m_GammaEvt) detID |= 0x80;
      addWord(detID,8,m_evtBufferCAL[tower],FES_CAL_MRLENG,&HeaderByteOff, &HeaderBitOff);

/* Now we have to add the zero hit transition record */
/* Length of zero hit record is 213 bytes*/
      int blkTran = 213;
// Start record header...Header has 5 bytes   
      recHead = ((blkTran)<<16) | (fesVersion<<8) | FES_RTYP_CAL_TRAN;
      addWord(recHead,32,m_evtBufferCAL[tower],FES_CAL_MRLENG,&byteOff, &bitOff);
      addWord(detID,8,m_evtBufferCAL[tower],FES_CAL_MRLENG,&byteOff, &bitOff);

      int time   = 40;      
      if(fesVersion==2) addWord(time,32,m_evtBufferCAL[tower],FES_CAL_MRLENG,&byteOff, &bitOff);
      
/* Now Add the transition vectors 4-32bit words = 16 Bytes */
      int datpt  = 0x1;
      int treqEmpty   = 0x0;
      unsigned int tVec = ((time & 0x3fff)<<16) | ((datpt & 0x7)<<8) | (treqEmpty & 0xff);   
      addWord(tVec,32,m_evtBufferCAL[tower],FES_CAL_MRLENG,&byteOff, &bitOff);
      addWord(tVec,32,m_evtBufferCAL[tower],FES_CAL_MRLENG,&byteOff, &bitOff);
      addWord(tVec,32,m_evtBufferCAL[tower],FES_CAL_MRLENG,&byteOff, &bitOff);
      addWord(tVec,32,m_evtBufferCAL[tower],FES_CAL_MRLENG,&byteOff, &bitOff);   

/* Now the blank cable summary words (6 32 bit words = 24 bytes) times four cables (96 bytes)*/
      for(int i=0; i<4; i++) {
         addWord(0,32,m_evtBufferCAL[tower],FES_CAL_MRLENG,&byteOff, &bitOff);
         addWord(0,32,m_evtBufferCAL[tower],FES_CAL_MRLENG,&byteOff, &bitOff);
         addWord(0,32,m_evtBufferCAL[tower],FES_CAL_MRLENG,&byteOff, &bitOff);
         addWord(0,32,m_evtBufferCAL[tower],FES_CAL_MRLENG,&byteOff, &bitOff);
         addWord(0,32,m_evtBufferCAL[tower],FES_CAL_MRLENG,&byteOff, &bitOff);
         addWord(0,32,m_evtBufferCAL[tower],FES_CAL_MRLENG,&byteOff, &bitOff);   

/* Pad out the data...to give an integer number of 8-48 bit words...we need 4x48 bit word = 8x24 bit words (24 bytes/cable) (96 bytes total)*/
         for(int i=0; i<24; i++) addWord(0,24,m_evtBufferCAL[tower],FES_CAL_MRLENG,&byteOff, &bitOff);
      }

/* Write out the buffer */
//     printData(beg,m_evtBufferCAL[tower],80);
//     writeCAL(tower,beg,nBytes+blkTran);      

// Save the information about the data to write after we know the time until the next event
// reset the buffer pointer to the beginning of the buffer.
     bLengthCAL[tower] = nBytes+blkTran;
     m_evtBufferCAL[tower] = beg;
     

   }// End Loop over towers 

   m_count++;
   if (m_count > 64) m_count=0;
      
   return 0;
}

/**
 *
 *  @fn   unsigned int FESOutput::completeCAL (int nDelta)
 *  @brief finishes write out of CAL FES data
 *  @param nDelta  number of clock ticks since last event
 *  @return     void
 *
 */
void FESOutput::completeCAL(int nDeltaTime){

// Loop over towers
  for(int tower=0; tower<16; tower++){

// If this is version 2, push in the delta time (LAT clock ticks to next event)
     if(fesVersion==2) {
         unsigned int timeWd = (1 << 31) | (nDeltaTime & 0x3fffffff);
//         printf("Event Time Word 0x%8.8x , EvtTimebyteOffCAL %i  EvtTimebitOff %i \n",timeWd,EvtTimebyteOffCAL,EvtTimebitOffCAL);
         int byteOff = EvtTimebyteOffCAL;
         int bitOff = EvtTimebitOffCAL;
         addWord(timeWd,32,m_evtBufferCAL[tower],FES_CAL_MRLENG,&byteOff,&bitOff);
     }   
     
// Write the file.
//     printf("Writing CAL Tower %i Buffer Length %i\n",tower,bLengthCAL[tower]);     
     writeCAL(tower,m_evtBufferCAL[tower],bLengthCAL[tower]);
     
// Keep the number of events written into this file. At the end of the job
// we will rewind the file and put the info in the header.
     m_CAL_evtCount[tower]++;     
  }

  return;
}

/**
 *
 *  @fn   unsigned int FESOutput::dumpACD (const EbfAcdData *acd)
 *  @brief Converts the detector data blocks into EBF format.
 *  @param tkr  The ACD data block
 *  @return     ???
 *
 */
unsigned int FESOutput::dumpACD(const EbfAcdData *acd, int nDeltaTime){


// Complete the time stamp on the header and write previous event to disk
   if(!firstEvt) completeACD(nDeltaTime);
      
   bool debug = false;
   if (debug) printf("FES dumpACD Routine\n");

// Associate the global cable number with corner and FES cable 
   int cornerCable[4][4]= { 1, 2,   3, -1,
                            4, 5,   6, -1,
                            0, 10, 11, -1,
                            7,  8,  9, -1 };

// Loop over corners/files
   for(int corner=0; corner<4; corner++) {
   
      if(debug) printf("Corner File %i \n",corner);
         
/*  The file header has already been written, the
    follow corresponds to the data record
    
      Record Header:  
      
      Detector(1byte)
      Cable list(4bytes)
          Up to 4 cables: Numbered 0-11          
                    
      Data(Variable Length)
          Transition Vector, cable 0  (6 bytes: 48 bits )
             bits 0-17  Vetoes Feed to GEM (Bit 0 => Chn 0 ; Bit 17 => Chn 17)
             bit    18  CNO Trigger
             bit 19-21  Control Advance of data memory?
             bit 22-31  Unused
             bit 32..45 Timer
             bit  46    Event Flag
             bit  47    Unused
          
          Transition Vector, cable 1   (6 bytes, if present)
          
          
          Transition Vector, cable 2   (6 bytes, if present)
          
          
          Transition Vector, cable 3    (6 bytes, if present)

 
          Hit Data, Cable 0  (variable length)
             Organized in 16 bit Words
             Wrd 0:   Chn 3 (4bits)  Chn 2 (4bits) Chn1 (4bits)  Chn 0 (4 bits)
             Wrd 1:   Chn 7          Chn 6         Chn5          Chn 4
             Wrd 2:   Chn 11         Chn 10        Chn 9         Chn 8
             Wrd 3:   Chn 15         Chn 14        Chn 13        Chn 10
             Wrd 4:   empty          empty         Chn 17        Chn 16
             
             For each Channel 4 bit nibble
                  bit 0 PHA word available [This is locally for the FES]
                  bit 1 A Hit is Available and above threshold...this is tied
                        to the Veto Bits for GLEAM.
                  bit 2 Accept Bit set for this channel
                  bit 3 Range Bit (0: Low Range; 1: High Range)
             
             Wrd 5:  ADC Count of first hit
             Wrd 6:  ADC Count of 2nd Hit
               .
               .
               .
             Wrd N:  ADC Count of Nth Hit (Max: 18)
             
             This data block must come in blocks of 8x32 Bits.  If the
             data is a maximum size, it consists of 5+18=23 16 bit words
             or 11.5 32 bit words...so this must be padded out to two
             blocks (8x32).  The minimum length is the 5 words for the 
             channel information.  This must be padded out to 1 8x32bit
             block. This padding is performed on a cable be cable basis.
 
          
          Hit Data, Cable 1  (variable length, if present)
 
          
          Hit Data, Cable 2  (variable length, if present)

          
          Hit Data, Cable 3  (variable length, if present)

*/

// Setup the buffer
      size_t buffSize = 2*FES_ACD_MRLENG;  
      memset(m_evtBufferACD[corner],0,buffSize);
      unsigned int *beg = m_evtBufferACD[corner];


// Add Blank word for Record header...must be filled at the end when we
//  know the length of the record.
     int byteOff = 0;
     int bitOff  = 0;
     int HeaderByteOff=byteOff;
     int HeaderBitOff =bitOff;
     addWord(0,32,m_evtBufferACD[corner],FES_ACD_MRLENG,&byteOff,&bitOff);
     
// Add byte for detector type
     int detID = FES_DETC_ACD | corner;
     if(m_GammaEvt) detID |= 0x80;
     addWord(detID,8,m_evtBufferACD[corner],FES_ACD_MRLENG,&byteOff,&bitOff);
     if(debug) printf("Detector ID Word 0x%8.8x \n",detID);

// Add 4 bytes for cable numbers
     int cableList = (cornerCable[corner][3]<<24)|(cornerCable[corner][2]<<16)|(cornerCable[corner][1]<<8)|cornerCable[corner][0];
     addWord(cableList,32,m_evtBufferACD[corner],FES_ACD_MRLENG,&byteOff,&bitOff);
     if(debug) printf("Cable List 0x%8.8x\n",cableList);


     EvtTimebyteOffACD =byteOff;
     EvtTimebitOffACD  =bitOff;
     if(fesVersion==2) addWord(0,32,m_evtBufferACD[corner],FES_ACD_MRLENG,&byteOff,&bitOff);
     
     int TVecByteOff = byteOff;
     int TVecBitOff  = bitOff;
     unsigned int TVec[4][2];

// Add blank words for cable transition vectors (6 bytes if present)
     for(int FESCable=0; FESCable<4; FESCable++){
        if(cornerCable[corner][FESCable]>=0) {
           addWord(0,32,m_evtBufferACD[corner],FES_ACD_MRLENG,&byteOff,&bitOff);
           addWord(0,16,m_evtBufferACD[corner],FES_ACD_MRLENG,&byteOff,&bitOff);
        }   
     }

// Loop over FES Cable
       for(int FESCable=0; FESCable<4; FESCable++){
          if(cornerCable[corner][FESCable]>=0){

// Get the EbfAcd Board Number
            int brd = cornerCable[corner][FESCable];

// Get the information from the FREE board
            unsigned int accepts = acd->BrdAccepts(brd);
            unsigned int vetoes  = acd->BrdVetoes(brd);       
//            printf("FES: Board %i Vetos 0x%8.8x  Acc  0x%8.8x\n",brd,vetoes,accepts);
            
// Begin Packing the Transition Vector words (Vetoes and CNO first)
// Must reverse the mapping of veto bits
//    bit 0 = Chn 0
//    bit 17 = Chn 17  (Ebf Quantity has this reversed
            TVec[FESCable][0] = 0;
            for(int chn = 0; chn<18; chn++){
               if((vetoes & (0x1 << (17-chn) ) ) > 0) TVec[FESCable][0] |= (0x1 << chn);
            }
            TVec[FESCable][0]|= ((acd->cno() >> cornerCable[corner][FESCable]) &0x1) <<18;

// Put in the timer value (0 for true event) and event flag
            int time =0;
            TVec[FESCable][1] = (1 << 14) | (time & 0x3fff);  
        
// Cycle through the accepts and count the number of hits packing
            int NHit =0;
            unsigned short adcVal[18];
            for(int chnl=0; chnl<acd->NumChannelsPerBoard; chnl++) {
            
               unsigned int wrd = 0;        
               if( (accepts & (0x1 << (17-chnl) ) ) > 0)  {
                  adcVal[NHit] = acd->BrdAdcs(brd,chnl);
                  unsigned int rng = acd->BrdAdcsRng(brd,chnl);   
                  wrd = ((rng&0x1)<<3)|(1<<2)|(0x1);
// If VETO is set for this channel set "HIT/VETO" bit
                  if ( (vetoes & (0x1 << (17-chnl) ) ) > 0) wrd |= 0x2;
                  if(debug) printf("NHit: %i Brd: %i Veto 0x%8.8x Acc 0x%8.8x Rng 0x%8.8x Chn 0x%8.8x Adc 0x%8.8x Wrd 0x%8.8x\n",NHit,brd,vetoes,accepts,rng,chnl,adcVal[NHit],wrd);
                  NHit++;
               }

// Add the nibble for the cable                 
               addWord(wrd,4,m_evtBufferACD[corner],FES_ACD_MRLENG,&byteOff,&bitOff);
               if(debug) printf("Nibble for Channel %i: 0x%8.8x\n",chnl,wrd);             

            } //loop over channels

// Pad out the channel information (empty byte)
            addWord(0,8,m_evtBufferACD[corner],FES_ACD_MRLENG,&byteOff,&bitOff);             


// Now fill the PHA hit values;
            for(int nhit=0; nhit<NHit; nhit++){
               addWord(adcVal[nhit],16,m_evtBufferACD[corner],FES_ACD_MRLENG,&byteOff,&bitOff);
               if(debug) printf("PHA Value #%i Stored: 0x%5.5x\n",nhit,adcVal[nhit]);
            }

// Pad out the Cable Data
            int n16BitWds = 5+NHit;
            int nPad = 16 - n16BitWds%16;
            for(int np=0; np<nPad; np++) addWord(0,16,m_evtBufferACD[corner],FES_ACD_MRLENG,&byteOff,&bitOff);            
            if(debug) printf("Padding Cable %i from %i 16 bit wds to %i 16 bit words\n",FESCable,n16BitWds,n16BitWds+nPad);

// Set the bits for the data memory advance in int Transition Vector
            int datp = n16BitWds <=16 ? 0x1 : 0x2;
            TVec[FESCable][0] |= (datp<<19);
            
          } //if we have an active cable
       }// loop over FESCable

// Put the transition vectors into buffer(6 bytes if present)...add upper bits first
       for(int FESCable=0; FESCable<4; FESCable++){
          if(cornerCable[corner][FESCable]>=0) {
            addWord(TVec[FESCable][0],32,m_evtBufferACD[corner],FES_ACD_MRLENG,&TVecByteOff,&TVecBitOff);
            addWord(TVec[FESCable][1],16,m_evtBufferACD[corner],FES_ACD_MRLENG,&TVecByteOff,&TVecBitOff);
            if(debug) printf("Cable %i TVec(48-33) 0x%8.8x TVec(0-32) 0x%8.8x\n",FESCable,TVec[FESCable][1],TVec[FESCable][0]);
          }   
       }
     
// Complete the record header.
      int nBytes = byteOff*4 + bitOff/8;  
      int recHead = ((nBytes)<<16) | (fesVersion<<8) | FES_RTYP_ACD_TRAN;
      addWord(recHead,32,m_evtBufferACD[corner],FES_ACD_MRLENG,&HeaderByteOff,&HeaderBitOff);
      if(debug) printf("Record Header (0x%8.8x Bytes): 0x%8.8x\n",nBytes,recHead);

// Add a blank Event Record

// (First Record header (5 bytes), then Detector ID)
       int HeadByteOff = byteOff;
       int HeadBitOff  = bitOff;
       addWord(0,32,m_evtBufferACD[corner],FES_ACD_MRLENG,&byteOff, &bitOff);
       addWord(detID,8,m_evtBufferACD[corner],FES_ACD_MRLENG,&byteOff, &bitOff);

// Cable List (4 bytes)
       addWord(cableList,32,m_evtBufferACD[corner],FES_ACD_MRLENG,&byteOff,&bitOff);

       int time = 40;
       if(fesVersion==2) addWord(time,32,m_evtBufferACD[corner],FES_ACD_MRLENG,&byteOff,&bitOff);
       
// Transition Vector (6 bytes each)  
       int nTran=0;                
       int unsigned TVecEmpty[2];
       TVecEmpty[0]=(0x1 << 19);
       TVecEmpty[1]=(time & 0x3fff);
       for(int FESCable=0; FESCable<4; FESCable++){
          if(cornerCable[corner][FESCable]>=0) {
            addWord(TVecEmpty[0],32,m_evtBufferACD[corner],FES_ACD_MRLENG,&byteOff,&bitOff);
            addWord(TVecEmpty[1],16,m_evtBufferACD[corner],FES_ACD_MRLENG,&byteOff,&bitOff);
            nTran++;
          }   
       }

// Add empty Data this is an 8x32 bit block of zeros, it represents
// the empty channel information and padding. (32 bytes total)
       for(int FESCable=0; FESCable<4; FESCable++){
          if(cornerCable[corner][FESCable]>=0) {
            for(int np=0; np<8; np++) addWord(0,8,m_evtBufferACD[corner],FES_ACD_MRLENG,&byteOff,&bitOff);
          }
       }
          
// Complete the record Header
       int blkData = 5+4+6*nTran+32*nTran;
       recHead = ((blkData)<<16) | (fesVersion<<8) | FES_RTYP_ACD_TRAN;
       addWord(recHead,32,m_evtBufferACD[corner],FES_ACD_MRLENG,&HeadByteOff, &HeadBitOff);


/* Write out the buffer */
//       printData(beg,m_evtBufferACD[corner],80);
//       writeACD(corner,beg,nBytes+blkData);     

// Save the information about the data to write after we know the time until the next event
// reset the buffer pointer to the beginning of the buffer.
     bLengthACD[corner] = nBytes+blkData;
     m_evtBufferACD[corner] = beg;
                                   
   } //Loop over corners (ACD FES Files)


return 0;
}

/**
 *
 *  @fn   unsigned int FESOutput::completeACD (int nDelta)
 *  @brief finishes write out of ACD FES data
 *  @param nDelta  number of clock ticks since last event
 *  @return     void
 *
 */
void FESOutput::completeACD(int nDeltaTime){

// Loop over towers
  for(int cnr=0; cnr<4; cnr++){

// If this is version 2, push in the delta time (LAT clock ticks to next event)
     if(fesVersion==2) {
         unsigned int timeWd = (1 << 31) | (nDeltaTime & 0x3fffffff);
//         printf("Event Time Word 0x%8.8x , EvtTimebyteOffCAL %i  EvtTimebitOff %i \n",timeWd,EvtTimebyteOffCAL,EvtTimebitOffCAL);
         int byteOff = EvtTimebyteOffACD;
         int bitOff = EvtTimebitOffACD;
         addWord(timeWd,32,m_evtBufferACD[cnr],FES_ACD_MRLENG,&byteOff,&bitOff);
     } 
     
// Write the file       
     writeACD(cnr,m_evtBufferACD[cnr],bLengthACD[cnr]);

// Keep the number of events written into this file. At the end of the job
// we will rewind the file and put the info in the header.
     m_ACD_evtCount[cnr]++;     
     
  }

  return;
}

/**************************************************************************
 **
 **  Add variable-length word of data to general data buffer
 **
 **  Description:
 **  -----------
 **
 **  Adds a specified number of bits of data to the next available bit
 **  location in a buffer.  The buffer is filled in a little-endian
 **  manner, so that at most 8 bits at a time are added to it.
 **
 **  Arguments:
 **  ---------
 **
 **    word     An integer containing the data to be added to the buffer.
 **
 **    wSize    An integer specifying the number of low-order bits of
 **             word to be added to the buffer.
 **
 **    buff     The address of the buffer to receive the data bits.
 **
 **    lBuff    The number of bytes in the buffer.
 **
 **    byteOff  The address of the byte offset of the current buffer
 **             location.  Updated by this routine call.
 **
 **    bitOff   The address of the bit offset within the current buffer
 **             byte.  Updated by this routine call.
 **
 **  Return Values:
 **  -------------
 **
 **    OK:     Success - the buffer can contain the data
 **    ERROR:  Failure - the buffer is too short for the data
 **
 **************************************************************************
 */
static int addWord(int word, int wSize, unsigned int *buff, int lBuff,
                   int *byteOff, int *bitOff)
{
    int size, bOff, nBits, cWord, mask;
    unsigned int *buffp;
    unsigned int value;

    buffp = buff + (byteOff ? *byteOff : 0);
    bOff = bitOff ? *bitOff : 0;
    cWord = word;
//    printf("ADD WORD: buff 0x%8.8x  byteOff %i   bitOff %i  data Value: 0x%8.8x\n",buff,*byteOff,*bitOff,word);

    for (nBits = wSize; nBits > 0; nBits -= size) {
        size = 32 - bOff;
        if (size > nBits) size = nBits;
        if (lBuff && buffp - buff >= lBuff) return 1;
//        if (bOff == 0) {
//            if (lBuff && buffp - buff >= lBuff) return 1;
//            value = 0;
//        }
//        else 
        value = *buffp;
        mask = size<32 ? (1 << size) - 1 : 0xffffffff;
   
        *buffp = value | ((cWord & (mask)) << bOff);
//        printf("bOff %i size %i mask 0x%8.8x value 0x%8.8x cWord 0x%8.8x buffp 0x%8.8x  Data 0x%8.8x\n",bOff,size,mask,value,cWord,buffp,*buffp);
        cWord >>= size;
        bOff += size;
        if (bOff == 32) {
            buffp++;
            bOff = 0;
        }
    }

    if (byteOff) *byteOff = buffp - buff;
    if (bitOff) *bitOff = bOff;

    return 0;
}

/**************************************************************************
 **
 **  Get variable-length word of data from general data buffer
 **
 **  Description:
 **  -----------
 **
 **  Gets a specified number of bits of data from the next available bit
 **  location in a buffer.  The buffer is assumed to be filled in a
 **  little-endian manner, so that at most 8 bits at a time are obtained
 **  from it.
 **
 **  Arguments:
 **  ---------
 **
 **    word     The address of an integer to contain the retrieved data.
 **
 **    wSize    An integer specifying the number of low-order bits of
 **             word to be obtained from the buffer.
 **
 **    buff     The address of the buffer containing the data.
 **
 **    lBuff    The number of bytes in the buffer.
 **
 **    byteOff  The address of the byte offset of the current buffer
 **             location.  Updated by this routine call.
 **
 **    bitOff   The address of the bit offset within the current buffer
 **             byte.  Updated by this routine call.
 **
 **  Return Values:
 **  -------------
 **
 **    OK:     Success - the data was retrieved
 **    ERROR:  Failure - the end of the buffer was reached
 **
 **************************************************************************
 */
static int getWord(int *word, int wSize, unsigned int *buff, int lBuff,
                   int *byteOff, int *bitOff)
{
    int size, bOff, nBits, value, vOff, mask;
    unsigned int *buffp;

    buffp = buff + (byteOff ? *byteOff : 0);
    bOff = bitOff ? *bitOff : 0;
    value = 0;
    vOff = 0;

    for (nBits = wSize; nBits > 0; nBits -= size) {
        size = 32 - bOff;
        if (size > nBits) size = nBits;
        if (bOff == 0 && lBuff && buffp - buff >= lBuff) return 1;
        mask = size<32 ? (1 << size) - 1 : 0xffffffff;
        value |= ((*buffp >> bOff) & (mask)) << vOff;
        vOff += size;
        bOff += size;
        if (bOff == 32) {
            buffp++;
            bOff = 0;
        }
    }

    *word = value;
    if (byteOff) *byteOff = buffp - buff;
    if (bitOff) *bitOff = bOff;

    return 0;
}

/**
 *
 *  @fn     unsigned int FESOutput::write ()
 *  @brief  Writes the current event to the output stream.
 *  @return The number of 32 bit words written
 *
 */


int * FESOutput::writeTKR (int tower, unsigned int *buff, int nBytes)
{
   return (int *)fwrite (buff, 1, nBytes, m_fpTKR[tower]);
}

int * FESOutput::writeCAL (int tower, unsigned int *buff, int nBytes)
{
   return (int *)fwrite (buff, 1, nBytes, m_fpCAL[tower]);
}

int * FESOutput::writeACD (int corner, unsigned int *buff, int nBytes)
{
   return (int *)fwrite (buff, 1, nBytes, m_fpACD[corner]);
}

/**
 *
 * @fn     int FESOutput::close ()
 * @brief  Closes the file associated with this output stream.
 * @return Status
 *
 */
unsigned int FESOutput::close ()
{
   int status;

// Make sure last event is written out
   if(!firstEvt){
      completeCAL(0);
      completeTKR(0);
      completeACD(0);
   }

//   printf("Closing the FES Files\n");
   /* If not already closed */
   for(int tower=0; tower<16; tower++){
     if (m_fpTKR[tower]) {

// For each file, rewind and look at top of file
           rewind(m_fpTKR[tower]);
           char wrds[10];
              
// Read to point where event/transition count should be stored
          unsigned int val = fread (wrds, sizeof (*wrds), 10, m_fpTKR[tower]);

//          printf("Event count %i \n",m_CAL_evtCount[tower]);
          unsigned int ntran = 2*m_TKR_evtCount[tower];
          for(int nbyte=0; nbyte<4; nbyte++) fputc(((ntran>>(8*nbyte))&0xff),m_fpTKR[tower]);
          for(int nbyte=0; nbyte<4; nbyte++) fputc(((m_TKR_evtCount[tower]>>(8*nbyte))&0xff),m_fpTKR[tower]);

// Close the file
          status = fclose (m_fpTKR[tower]);
          m_fpTKR[tower] = 0;
     } 
   } 

   /* If not already closed */
   for(int tower=0; tower<16; tower++){
     if (m_fpCAL[tower]) {
     
// For each file, rewind and look at top of file
           rewind(m_fpCAL[tower]);
           char wrds[10];
              
// Read to point where event/transition count should be stored
          unsigned int val = fread (wrds, sizeof (*wrds), 10, m_fpCAL[tower]);

//          printf("Event count %i \n",m_CAL_evtCount[tower]);
          unsigned int ntran = 2*m_CAL_evtCount[tower];
          for(int nbyte=0; nbyte<4; nbyte++) fputc(((ntran>>(8*nbyte))&0xff),m_fpCAL[tower]);
          for(int nbyte=0; nbyte<4; nbyte++) fputc(((m_CAL_evtCount[tower]>>(8*nbyte))&0xff),m_fpCAL[tower]);

// close the file               
          status = fclose (m_fpCAL[tower]);
          m_fpCAL[tower] = 0;
     } 
   }
   /* If not already closed */
   for(int corner=0; corner<4; corner++){
       if (m_fpACD[corner]) {
       
// For each file, rewind and look at top of file
           rewind(m_fpACD[corner]);
           char wrds[10];
              
// Read to point where event/transition count should be stored
          unsigned int val = fread (wrds, sizeof (*wrds), 10, m_fpACD[corner]);

//          printf("Event count %i \n",m_CAL_evtCount[tower]);
          unsigned int ntran = 2*m_ACD_evtCount[corner];
          for(int nbyte=0; nbyte<4; nbyte++) fputc(((ntran>>(8*nbyte))&0xff),m_fpACD[corner]);
          for(int nbyte=0; nbyte<4; nbyte++) fputc(((m_ACD_evtCount[corner]>>(8*nbyte))&0xff),m_fpACD[corner]);
       
// Close the file       
            status = fclose (m_fpACD[corner]);
            m_fpACD[corner] = 0;
       }       
   }
/* FEATURE: Only status on last closing...need to fix */   
   return status;
}






/**
 *  
 *  @fn    FESOutput::~FESOutput ()
 *  @brief Closes the file associated with this stream (if it hasn't 
 *         already been closed) and frees the internal memory it 
 *         has allocated.
 */
FESOutput::~FESOutput ()
{
   /* If have anything left in the event buffer, write it out */
  //if (m_curEvtSize) write ();
 
 
   /*
    | Must close the file if it is still opened and free the memory
    | associated with the event buffers, numbers and lengths. All
    | three of these blocks were allocated contigiously so just need
    | to free the memory associated with the top of this block. This
    | is the event buffer.
   */
//   printf("In the destructor of FESOutput\n");
//   close ();
/*
   free (m_evtBufferTKR);
   free (m_evtBufferCAL);
   free (m_evtBufferTKR);
*/   
}





/**
 *
 *  @fn    void FESOutput::print () const
 *  @brief Prints the event(s) currently in the output buffer to 
 *         standard output stream. 
 *
 */
void FESOutput::print () const
{

    
    return;
}

static void printData (const unsigned int *beg, 
                       const unsigned int *end,
                       int         rightMargin)
/*
   DESCRIPTION
   -----------
   Prints the specified binary data in the form of a hex dump. This is
   intentionally crude, and is used primarily as a low level debugging
   aid.

   PARAMETERS
   ----------
          beg: Pointer to the beginning of the data to print.

          end: Pointer to the ending    of the data to print.

  rightMargin: Column number of the right margin or, alternatively,
               the maximum length of a line.

  RETURNS
  -------
  Nothing
*/
{
    int ncol = 0;
  
    for (const unsigned int *p = beg; p < end; p++)
    {
        /* If no room for the next value, emit a line feed */
        if (ncol > rightMargin) 
        {
            std::printf ("\n");
            ncol = 0;
        }
           
        ncol += std::printf (" %8.8x", *p);
    }

    std::printf ("\n");

    return;
}
