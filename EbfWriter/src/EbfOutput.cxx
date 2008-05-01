#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>

#include "EbfAcdData.h"
#include "EbfCalData.h"
#include "EbfTkrData.h"
#include "EbfGemData.h"
#include "EbfOutput.h"
#include "EbfContrib.h"



#undef DEBUG
#ifdef DEBUG
#define _DBG(statement) statement
#else
#define _DBG(statement)
#endif

/*
 *
 * The data from all contributors is packaged by LCB. LCB, itself
 * always transports its data in an integer number of LCB packets.
 * If the contributor data is not an integer number of LCB packets,
 * LCB pads the data until the data is an integer number of LCB
 * packets. LCB_PACKET_SIZE defines the size of an LCB packet.
 *
 */
#define LCB_PACKET_SIZE 16 /* LCB PACKET SIZE is 16 bytes */

/*
 *
 * The data from the event must be packed in 128 byte packages
 * when going through the GASU/EBM. If the event is short of
 * this then must pad out the end.
 *
 */
#define EVT_PACKET_SIZE 128 /* EVT PACKET SIZE is 128 bytes (GASU/EBM) */

/* Define a dummy EPU to act as the destination */
#define EPU 0x20

/* Define the size in 32 bit words of the event header (before contributions) */
#define EVT_HEADER_SIZE 8 
#define LDF_HEADER_SIZE 4

/* Size of the circular buff in the LCB in Bytes 512*1024 = 524288 = 2^19 */
/* Circular Buffer offset is in units of 32 bit words (4 bytes), so divide the */
/* size by 4:  2^19 --> 2^17 = 131072 */
#define CIRC_BUFFER_SIZE 131072

/* Maximum size of the EBF contribution 4096 - 128 + 32 = 4000 bytes */
#define EBF_MAX_PACKET_SIZE 4000
// For Testing
//#define EBF_MAX_PACKET_SIZE 160    /* 1*128 + 32 */
//#define EBF_MAX_PACKET_SIZE 288  /* 2*128 + 32 */
//#define EBF_MAX_PACKET_SIZE 800  /* 6*128 + 32 */
//#define EBF_MAX_PACKET_SIZE 1568 /* 12*128+32 */


/* This is the standard header placed on all LCB data packets. */
struct LcbHeader
{
    unsigned int hdr_len;
    unsigned int summary;
};


    
static inline unsigned int   *_advance (const unsigned int *ptr,
                                        int              nbytes);
static inline unsigned int  headerWord (unsigned int sourceId,
                                        unsigned int destinationId);
static inline unsigned int summaryWord (unsigned int sequence,bool FourRange);
static inline unsigned int         pad (const unsigned int *beg,
                                              unsigned int *end);

static inline unsigned int         padEvt(const unsigned int *beg,
                                                unsigned int *end);
                                              
static inline unsigned int   *complete (unsigned int *beg,
                                        unsigned int *end,
                                        unsigned int contributorId, 
                                        unsigned int sequence,
                                        unsigned int length);

static inline void                swap (unsigned int *wrds,
                                        int          nwrds);
static inline unsigned int      swrite (FILE         *fp,
                                        unsigned int *wrds,
                                        int          nwrds);


static void        printOneContributor (const unsigned int *beg,
                                        const unsigned int *end,
                                        int         rightMargin);
static int                 printOneEvt (const unsigned int *evt);
static void                  printData (const unsigned int *beg, 
                                        const unsigned int *end,
                                        int         rightMargin);




static unsigned int *_advance (const unsigned int *ptr,
                              int               nbytes)
/*
   DESCRIPTION
   -----------
   Advances the integer pointer by the specified number of bytes.
   This routine exists for the sole purpose of hiding the casts
   that need to be made.

   PARAMETERS
   ----------
          ptr: The pointer to advance.

       nbytes: The number of bytes to advance the pointer.
   
   RETURNS
   -------
   The advanced pointer.
 */
{
  return (unsigned int *)((unsigned char *)ptr + nbytes);
}



static unsigned int pad (const unsigned int *beg,
                               unsigned int *end)
/*
   DESCRIPTION
   -----------
   Pads the LCB packet that extends from 'beg' to 'end' to 
   end on an integer number of LCB packets. This usually involves
   adding a few 0-filled bytes to the packet.

   PARAMETERS
   ----------
          beg: Pointer to the beginning of the packet

          end: Pointer to the current end of the packet

   RETURNS
   -------
   The length of the padded packet in bytes.
*/
{
   unsigned int nbytes = (unsigned char *)end - (unsigned char *)beg;
   int            npad = nbytes % LCB_PACKET_SIZE;

   /* If not an integral number of packets, round up */
   if (npad)
   {
       unsigned char *dst = (unsigned char *)end;
       npad    = LCB_PACKET_SIZE - npad;
       nbytes += npad;
       while (--npad >= 0) *dst++ = 0x00;
   }
   
   /* Return the size in bytes */ 
   return nbytes;
}

static unsigned int padEvt (const unsigned int *beg,
                               unsigned int *end)
/*
   DESCRIPTION
   -----------
   Pads the event that extends from 'beg' to 'end' to 
   end on an integer number of EVT packets. This usually involves
   adding a few 0-filled bytes to the packet.

   PARAMETERS
   ----------
          beg: Pointer to the beginning of the packet

          end: Pointer to the current end of the packet

   RETURNS
   -------
   The length of the padded packet in bytes.
*/
{
   unsigned int nbytes = (unsigned char *)end - (unsigned char *)beg;
   int            npad = ( (nbytes-4) % EVT_PACKET_SIZE);

   /* If not an integral number of packets, round up */
//   printf("Number of Bytes %i Remainder %i\n",nbytes,npad);
   if (npad)
   {
       unsigned char *dst = (unsigned char *)end;
       npad    = EVT_PACKET_SIZE - npad;
       nbytes += npad;
//       printf("npad bytes %i nbytes %i\n",npad,nbytes);
       while (--npad >= 0) *dst++ = 0x00;
   }
   
   /* Return the size in bytes */ 
   return nbytes;
}


static unsigned int headerWord (unsigned int sourceId,
                                unsigned int destinationId)
/*
   DESCRIPTION
   -----------
   Composes the header word from the LCB source and destination ids.
   This routine also adds in the correct parity for the composed
   header word.

   PARAMETERS
   ----------
     sourceId: A small integer giving the identity of the source,
               that is, whose data is this.

destinationId: A small integer giving the identity of the destination,
               that is, which EPU is this destined to.
     

   RETURNS
   -------
   The compose header word.
*/
{
   unsigned int parity;

   sourceId      &= 0x3f;
   destinationId &= 0x3f;
   unsigned int LATp = 1; // LATp Protocol.

   unsigned int LATpWrd = (1 << 15) | (destinationId << 9) | (LATp & 0x3) << 7 | (sourceId << 1);

   /* Compute the parity over the full LATp word */
   parity       =  LATpWrd ^ (LATpWrd>> 8);
   parity       =  parity  ^ (parity >> 4); 
   parity       =  parity  ^ (parity >> 2);
   parity       = (parity  ^ (parity >> 1)) & 1;
  
   /* The bit at 15 is the response bit which indicates that additional
   / contributions follow.  We set this for all and then zero the bit in
   / the last contribution when we know what it is.
   */
//   printf("sourceID 0x%8.8x destinationId 0x%8.8x parity %i\n",sourceId,destinationId,parity);
   return LATpWrd | (parity ^ 1);
}




static unsigned int summaryWord (unsigned int sequence, bool FourRange)
/*
  DESCRIPTION
  -----------
  Composes the summary word. Currently only those fields dealing with
  the event sequence are filled in. If and when the other fields, such
  as the 4-range readout and the threshold suppress are necessary, they
  should be added to this routine.

  For the most part, this word is exactly the trigger word.

  PARAMETERS
  ----------
    sequence: The event sequence number.

  RETURNS
  -------
  The summary word for this event.
*/
{
   
  int summary = 0;

// For now setting these by hand but should pull them from the MC Event
  int startBit  = 1;
  int calStrobe = 0;
  int trigAck   = 1;
  int fourRange = (FourRange) ? 1:0;
  int zeroSup   = 1;
  int marker    = 0;
  int ErrCont   = 0;
  int DiagCont  = 0;
  int ParErr    = 0;
  
// Set the bits in the summary word
  summary |= (startBit & 0x1) << 31;         //Set Start bit
  summary |= (calStrobe & 0x1)<< 30;         //Cal Strobe
  summary |= (sequence & 0x3) << 28;         //Lower two bits of event number
  summary |= (trigAck & 0x1)  << 27;         //Trigger Acknowlege
  summary |= (fourRange & 0x1)<< 26;         //Four Range Readout
  summary |= (zeroSup & 0x1)  << 25;         //Zero Suppress
  summary |= (marker & 0x7)   << 22;         // Marker Bits
  summary |= (ErrCont & 0x1)  << 21;         // Error Contribution Present
  summary |= (DiagCont & 0x1) << 20;         // Diagnostic Contribution Present
  summary |= ((sequence >> 2) & 0x7fff) <<1; // High Evt Count
  summary |= (ParErr & 0x1);                 // Parity Error
  
  return  summary;
}




static unsigned int *complete (unsigned int *beg, 
                               unsigned int *end,
                               unsigned int sourceId,
                               unsigned int destinationId,
                               unsigned int summary)
/*
   DESCRIPTION
   -----------
   Completes the LCB packet defined by 'beg' and 'end' pointers. Completing
   the packet involves two operation. The first is padding the packet 
   out to an integer number of LCB packets. The second is filling in the
   LCB header.

   PARAMETERS
   ----------
          beg: A pointer to the beginning of the packet. This will
               be where the LCB header is composed

          end: A pointer to the ending of the packet. This address
               will be extended, if necessary, so that the packet
               is an even number of bytes.

     sourceId: A small integer giving the identity of the source,
               that is, whose data is this.

destinationId: A small integer giving the identity of the destination,
               that is, which EPU is this destined to.

      summary: The event summary word. For the most part, this is
               identically the trigger message.

  RETURNS
  -------
  A pointer to where the next LCB packet will begin.
*/
{
   struct LcbHeader *hdr = (struct LcbHeader *)beg;
   int            length = pad (beg, end);

   /* Final Length Reported is in 128 bit (16 byte) words */
   int lcbWords = length/16;
   if(lcbWords > 0x255) printf("lcbWords > 255 on Contributor %i sourceId\n",sourceId);

   int err = 0; /* Error word...no meaning in GLEAM...hardware bits */
   int seq = 0; /* Bit field for truncated events...set to 0 as initial packet*/
   int hdrWord = headerWord (sourceId, destinationId);
   hdr->hdr_len = ( hdrWord << 16) | (err & 0x7) << 13 | (seq & 0x1f)<<8 | (lcbWords & 0xff);
   hdr->summary = summary;

   /*
   printf ("Filling in the header at %p\n"
           "hdr->hdr_len = %8.8x\n"
           "hdr->summary = %8.8x\n",
           hdr,
           hdr->hdr_len,
           hdr->header,
           hdr->summary);
   */

   return _advance (beg, length);
}





static void swap (unsigned int *wrds, int nwrds)
/*
   DESCRIPTION
   -----------
   Performs a byte-swap on the specified number of 32-bit words. This
   is an El-Cheapo implementation, it should not be used for heavy
   duty byte-swapping.

   PARAMETERS
   ----------
         wrds: The 32-bit words to byte-swap

        nwrds: The number of 32-bit words to byte-swap. 

   RETURNS
   -------
   Nothing
*/
{
   while (--nwrds >= 0)
   {
     unsigned int tmp;
     
     tmp = *wrds;
     tmp = ((tmp & 0x000000ff) << 24) |
       ((tmp & 0x0000ff00) <<  8) |
       ((tmp & 0x00ff0000) >>  8) |
       ((tmp & 0xff000000) >> 24); 
     *wrds++ = tmp;
   }
   
   return;
}






static unsigned int swrite (FILE           *fp,
                            unsigned int *wrds,
                            int          nwrds)
/*
    DESCRIPTION
    -----------
    Writes the specified number of 32-bit words to the specified output
    stream. If necessary, the 32-bit words are byte-swapped into 
    big-endian format. Note, if the buffer is byte-swapped, it is done
    in place.

    PARAMETERS
    ----------
            fp: The file handle of the output stream.

          wrds: The 32-bit words to write out.

         nwrds: The number of 32-bit words to write out.

   RETURNS
   -------
   The return value of the output routine.
*/
{
   /*
    |  !!! KLUDGE !!!
    |  --------------
    |  Need to know whether the machine executing this code is a big
    |  or little endian machine. I've checked around and found no one
    |  who can tell me of a compiler defined symbol containing this
    |  tidbit of information. Currently, I've only every run GLEAM on
    |  Intel processors, so I've hardwired this to be little-endian.
   */
    swap (wrds, nwrds);
    if(fp) return fwrite (wrds, sizeof (*wrds), nwrds, fp);
    return 0;
}






/**
 *
 *  @fn      EbfOutput::open (const char    *fileName,
 *                            unsigned int maxEvtSize)
 *  @brief Initializes the EbfOutput be opening the named file and allocating
 *         enough memory to hold at least 1 maximally sized event.
 *
 *  @param fileName   The name of the file to open.
 *  @param maxEvtSize The maximum size of a given event.
 *
 */
int EbfOutput::open (const char     *fileName,
                     unsigned int  maxEvtSize)
{
   unsigned char           *ptr;

//   printf("EbfOutput::open file %s\n",fileName);

   /* Initialize all the data members to something beign */   
   m_fp         = 0;
   m_maxEvtSize = maxEvtSize;
   m_curEvtSize = 0;
   m_totEvtSize = 0;
   m_numEvtsIn  = 0;
   m_numEvtsOut = 0;
   m_evtBuffer  = 0;
   m_circBuffOff= 0;

   
   /* Allocate a buffer big enough to hold a maximally sized event */
   ptr = (unsigned char *)malloc (maxEvtSize);
   if (ptr == 0) { return -1;  }
   m_evtBuffer = (unsigned int *)ptr;
   
   /* Open the output file */
   char file[120];
   if(LdfFormat){sprintf(file,"%s.ldf",fileName);}else{sprintf(file,"%s.ebf",fileName);}

   if( ::strcmp(fileName, "")!=0) {
     m_fp = fopen (file, "wb");
     if (!m_fp)
     {
       /* Error in opening the output file.. */
         throw std::invalid_argument("EbfOutput::Error cannot open Ebf output file");
       free (ptr);
       return -2;
     }
    /* Report the name of the output file, maximum event size */  
#if 0// printf statements not acceptable in Gaudi algorithms
     printf ("EbfOutput::initialize: FileName: %s\n"
     "                 Max Event Size: %d\n",
     fileName,
     maxEvtSize);
#endif
   } else {
#if 0 // printf statements not acceptable in Gaudi algorithms
      /* Report the name of the output file, maximum event size */  
      printf ("EbfOutput::initialize: No Output file defined \n");
#endif
   }

  
   return 0;
}





/**
 *
 *  @fn   unsigned int EbfOutput::format (const EbfAcdData *acd,
 *                                        const EbfCalData *cal,
 *                                        const EbfTkrData *tkr,
 *                                        const EbfGemData *gem)
 *  @brief Converts the detector data blocks into EBF format.
 *  @param acd  The ACD data block.
 *  @param cal  The CAL data block, contains all 16 towers.
 *  @param tkr  The TKR data block, contains all 16 towers.
 *  @param gem  The GEM data block.
 *  @param mcInfo  Monte Carlo information to pack into header of ebf
 *  @param reGemTrig Require Gem Trigger for output
 *  @return     The number of 32 bit words needed to format the event
 *
 */
unsigned int EbfOutput::format (const EbfAcdData *acd,
                                const EbfCalData *cal,
                                const EbfTkrData *tkr,
                                const EbfGemData *gem,
                                unsigned int *mcInfo,
				bool reqGemTrig)
{
    int print = m_print;
    

   /*  
    | This convention matches the event display's numbering convention,
    | ie. the first event out is labeled #1. This number will be tacked
    | on to the end of the GEM data (hidden away).
    */
   unsigned int numEvtsIn = m_numEvtsIn + 1;
   
   
   /* Count the number of events seen */
   m_numEvtsIn = numEvtsIn;
 

   /* Check if this event is triggered */
   if (!gem->isTriggered () && reqGemTrig)
   {
//       printf (" WARNING: Event in = %8d not triggered\n", numEvtsIn);
       return 0;
   }

   
   unsigned int numEvtsOut = m_numEvtsOut;


   /* Count the number of events output */   
   m_numEvtsOut = numEvtsOut + 1;

   /* Determine if this is a 4 Range Readout */
   bool FourRange = cal->FourRangeReadOut();

   /* Event is triggered, format it to the output buffer */
   unsigned int      summary = summaryWord (numEvtsOut,FourRange); // Event sequence #
   /* Start of prefix held onto through out buffer build*/
   unsigned int      *evtBeg = m_evtBuffer;             

   /* Next three points are reset at the beginning of each contributor */
   /* Start of contributor */
   unsigned int         *beg = (LdfFormat) ? evtBeg + LDF_HEADER_SIZE : evtBeg + EVT_HEADER_SIZE;   // Contributor begin
   /* Start of data after event builder word and event summary word (all part of LcbHeader Structure) */
   unsigned int               *dst;
   unsigned int *LastContrib = beg;   

   /* Pack up the GEM contribution */
   dst    = _advance (beg, sizeof (LcbHeader));
   dst    = gem->format (dst);
   dst    = complete (beg, dst, EbfContributorId::GEM, EPU, summary);
   //gem->print();
   /*
     if (print)
   {
     printf ("GEM data:\n"); 
     printData (beg, dst, 72);
     }
   */
   


   /* Fill in the Tower Contributions */
   for (unsigned int itower = 0; itower < EbfCalData::NumTowers; itower++)
   {
       beg = dst;
       LastContrib = beg;       
       dst = _advance (beg, sizeof (LcbHeader));
       dst = cal->format (dst, itower);
       dst = tkr->format (dst, itower);
       dst = complete (beg, dst, EbfContributorId::TWR+itower, EPU, summary);
       /*
	 if (print)
	 {
	 printf ("TWR %2d data:\n", EbfContributorId::TWR+itower);
	 printData (beg, dst, 72);
	 }
       */
       
   }


   /* Fill in the ACD Contribution */
   beg    = dst;
   LastContrib = beg;
   dst    = _advance (beg, sizeof (struct LcbHeader));
   dst    = acd->format (dst);
   dst    = complete (beg, dst, EbfContributorId::ACD, EPU, summary);
   /*
     if (print)
     {
     printf ("ACD data:\n"); 
     printData (beg, dst, 62);
     }
   */
   
   /* Now for the last contributor we need to zero out the response bit */
//   printf("Last Contrib: 0x%8.8x \n",*LastContrib);
   *LastContrib &= 0x7fffffff;
   /* Add flip the parity bit */
   *LastContrib ^= 0x10000;
//   printf("Last Contrib (Fixed): 0x%8.8x \n",*LastContrib);

   int contribLength = 0;
   if(!LdfFormat) {
   /* Fill the event header words (EVT_HEADER_SIZE 32 bit words) 
   /  
   /    Start by padding event out to integer number of EVT_PAD_LENGTH  
   /   Word 0: Application: McInfo (sequence number)
   /   Word 1:            : McInfo ( 0-15: Reserved  16-23: First AcdTile (0xaa) # 24:31: source ID)
   /   Word 2:            : McInfo ( 0-15: absRA (0xcccc) 16-31: absDec (0xdddd))
   /   Word 3:            : McInfo ( 0-15: relTheta 16-31: relPhi)
   /   Word 4:            : McInfo ( 0-15: MC Energy of Primary 16-31: Observed MC Energy (0x2222) )
   /   Word 5: FSW Reserved  (Currently stored as 0x1111111)
   /   Word 6: FSW Reserved  (Currently stored as 0x2222222)
   /   Word 7: Descriptor 
   /             Bit 0-16:  Circular Buffer Offset
   /             Bit 17-26: Size of Contributors in 32 bit words
   /             Bit 27-29: Transmit Status
   /             Bit 30-31: Receive Status
    */
//   unsigned int length = (unsigned char *)dst - (unsigned char *)evtBeg;
//   int contribLength = length - 4*(EVT_HEADER_SIZE-1);  
       *evtBeg++      = *mcInfo++;
       *evtBeg++      = *mcInfo++;
       *evtBeg++      = *mcInfo++;
       *evtBeg++      = *mcInfo++;
       *evtBeg++      = *mcInfo++;

       *evtBeg++      = 0x11111111;
       *evtBeg++      = 0x22222222;
//        contribLength = padEvt(evtBeg,dst)-4; /*Subtract off Descriptor, returned value is bytes*/         
        contribLength = (dst - evtBeg)*sizeof(dst) - 4; /*Subtract off Descriptor, returned value is bytes*/        
        m_evtDescriptor = evtBeg; /* Save position of event descriptor */
        m_evtEnd        = dst; /* Save position of the end of the event */
//        printf("Contribution Length %i (0x%8.8x)\n",contribLength/4,contribLength/4);
       if(m_circBuffOff > (m_circBuffOff & 0x1ffff) ) printf("WARNING: Circular Buffer Too Large\n");
       *evtBeg++      = ( (contribLength/4) << 17) | (m_circBuffOff & 0x1ffff);
        m_evtHead     = evtBeg;
  
//   printf("contribLength 0x%8.8x 0x%8.8x\n",evtBeg,contribLength); 
       m_curEvtSize  = contribLength+EVT_HEADER_SIZE*4;
//       m_totEvtSize += m_curEvtSize;
//       m_circBuffOff = m_totEvtSize%CIRC_BUFFER_SIZE; 

//      printf("Last Contrib (Last): 0x%8.8x \n",*LastContrib);

    } else {

// LDF Packaging of data:
/*
/       Word 0: LATdatagram Identity: 
/       Word 1: LATdatagram Length (includes header words)
/       Word 2: LATcontribution Identity
/       Word 3: LATcontribution Length (includes Contribution Id and Length Word)
*/
       int LATdatagramID = 0xF1010;
       int LD_Version    = 0x201;
       *evtBeg++ = (LD_Version)<<20 | (LATdatagramID & 0xfffff);    
       contribLength = padEvt(evtBeg,dst);
       *evtBeg++ = contribLength+(LDF_HEADER_SIZE-1)*4;
       
       int LATcontribID = 0xF0010;
       int LC_Version   = 0x104;
       *evtBeg++ = (LC_Version)<<20 | (LATcontribID & 0xfffff);
       *evtBeg++ = contribLength+(LDF_HEADER_SIZE-3)*4; 

// Note: the contribLength includes the last header word
       m_curEvtSize  = contribLength+(LDF_HEADER_SIZE-1)*4;
       m_totEvtSize += m_curEvtSize;
       
    }
   //printf("Total Event Size: 0x%8.8x\n",m_curEvtSize);
      
   /* Return the length */
   return m_curEvtSize;
}


/**
 *
 *  @fn   unsigned int EbfOutput::writeMC (unsigned int *mcInfo)
 *  @brief Writes a LATDatagram for the Monte Carlo Information.
 *  @param mcInfo  Pointer to Monte Carlo information 
 *  @param size    Size of MC block in bytes
 *
 */
void EbfOutput::writeMC (unsigned int *mcInfo, int size)
{

       unsigned int *evtBeg = m_evtBuffer;
       unsigned int *data   = evtBeg + LDF_HEADER_SIZE;


// Now add the monte carlo information
       for(int i=0; i<size/4; i++) *data++ = *mcInfo++;     
       m_curEvtSize = padEvt(evtBeg,data)-4; 
       
// LDF Packaging of data:
/*
/       Word 0: LATdatagram Identity: 
/       Word 1: LATdatagram Length (includes header words)
/       Word 2: LATcontribution Identity
/       Word 3: LATcontribution Length (includes Contribution Id and Length Word)
*/
       int LATdatagramID = 0xF1020;
       int LD_Version    = 0x201;
       *evtBeg++ = (LD_Version)<<20 | (LATdatagramID & 0xfffff);    
       *evtBeg++ = m_curEvtSize;   
       
       int LATcontribID = 0xF0020;
       int LC_Version   = 0x101;
       *evtBeg++ = (LC_Version)<<20 | (LATcontribID & 0xfffff);
       *evtBeg++ = m_curEvtSize - 2*4; /* subtract off datagram header */ 

         
// Write out the block
       unsigned int length;
       char *out = write(true,length,0);
        
   return;
}

/**
 *
 *  @fn     unsigned int EbfOutput::write ()
 *  @brief  Writes the current event to the output stream.
 *  @return The number of 32 bit words written
 *
 */
char * EbfOutput::write (bool writeEbf, unsigned int &length, unsigned int *TdsBuffer)
{
   unsigned int toWrite = m_curEvtSize;
   bool truncate = false;
   int seq = 0;
   
   int i;
   int offset = 0;
   length = 0;


// Check to see if the event is too big. If so, we must truncate
   do {

// This is only called on the subsequent passes through the
// do loop...it is for writing the remainder of the event       
       if(truncate) {

// Increment the sequence word
          seq++; 

// The part of the buffer which we will rewrite, needs to be
// byte swapped because it was byte swapped when writting out.
          if(writeEbf) swap(m_evtBuffer,9);
     
// Add prefix words and new LATp header word
//          printf("Modifying the LATpHeader\n");
          int LATpHeader = *m_evtHead;
//          printf("Existing LATpHeader 0x%8.8x \n",LATpHeader);
          LATpHeader = (LATpHeader & 0xffff0000) | (seq &0x1f) <<8 | 0x1;
//          printf("New LATpHeader 0x%8.8x \n",LATpHeader);
          unsigned int * evtData = m_evtHead;
          *evtData++ = LATpHeader;
// Add three zero words to fill out the 128 bit cell
          *evtData++ = 0;
          *evtData++ = 0;
          *evtData++ = 0;

// Move truncated data up in m_evtBuffer
// Assumes that EBF_MAX_PACKET_SIZE includes header
          unsigned int *evtRem = m_evtHead;
          for(int i=0; i<(EBF_MAX_PACKET_SIZE-EVT_HEADER_SIZE*4)/sizeof(evtRem); i++) {
//            printf("Pointer Position 0x%8.8x Value 0x%8.8x\n",evtRem,*evtRem);
            *evtRem++;
          }  
          
// Move data up
          int totalWd = (toWrite - EBF_MAX_PACKET_SIZE)/sizeof(*evtData);
//          printf("Shifting Data up totalWd %i\n",totalWd);
          for(int nwd=0; nwd<totalWd; nwd++) {
//             printf("nWd %i Old Word 0x%8.8x New Word 0x%8.8x Pointer 0x%8.8x\n",nwd,*evtData,*evtRem,evtRem);  
             *evtData++ = *evtRem++;
          }

// Keep track of the end of the event
          m_evtEnd = evtData;
//          printf("End of Data Pointer 0x%8.8x \n",m_evtEnd);
          
// Determine new length to write
//  We will rewrite the header words and we have three extra words
//  of the restart cell.
          int written = (EBF_MAX_PACKET_SIZE-EVT_HEADER_SIZE*4);
//          printf("Starting toWrite %i \n",toWrite);
          toWrite -= (EBF_MAX_PACKET_SIZE-EVT_HEADER_SIZE*4-4*4) ;
//          printf("New toWrite %i\n",toWrite);
       }
       
// Do we need to truncate the remainder of the event       
       truncate = (toWrite > EBF_MAX_PACKET_SIZE);
//       printf("Truncate %i Seq %i toWrite %i  EBF_MAX_PACKET_SIZE %i \n",truncate,seq,toWrite,EBF_MAX_PACKET_SIZE);
       unsigned int packetSize = toWrite;

// If we need to truncate this current packet, set the bits in the 
// event descriptor to indicate packet truncation.
       if(truncate & !LdfFormat) {
          
            unsigned int Descriptor = *m_evtDescriptor; 
//            printf("Event Desc 0x%8.8x...truncation\n",Descriptor);
// Make sure the RSTATUS and XSTATUS bits are zeroed, then set by hand.
            *m_evtDescriptor = (Descriptor & 0x0001ffff) | ((EBF_MAX_PACKET_SIZE-EVT_HEADER_SIZE*4)/4<<17) | (0x3 << 30);
//            printf("New Event Descriptor 0x%8.8x....truncation\n",*m_evtDescriptor);
            packetSize = EBF_MAX_PACKET_SIZE;
            

// Keep running count of circular buffer postion
            m_totEvtSize += packetSize;
// Total event size is in bytes. Circular Buffer Offset is in 32-bit words             
            m_circBuffOff = (m_totEvtSize/4)%CIRC_BUFFER_SIZE; 

       } else if (!truncate & !LdfFormat) {

// Make sure this last packet is padded to a 128 byte boundary. Adding padding if necessary
            unsigned int npad = EVT_PACKET_SIZE - (toWrite - EVT_HEADER_SIZE*4)%EVT_PACKET_SIZE;
            if(npad==EVT_PACKET_SIZE) npad=0;
//            printf("toWrite 0x%8.8x npad %i \n",toWrite,npad); 
//            unsigned char *dst = (unsigned char *)m_evtEnd;
            toWrite += npad;
//            while (--npad>=0) *dst++=0xFF;
            int nwd = npad/4;
//            printf("Number of words to pad %i \n",nwd);
            for(int i=0; i<nwd; i++) *m_evtEnd++ = 0x00000000;
//            printf("New toWrite 0x%8.8x \n",toWrite);
            packetSize = toWrite;

// Keep running count of circular buffer postion
            m_totEvtSize += packetSize; 
// Total event size is in bytes. Circular Buffer Offset is in 32-bit words             
            m_circBuffOff = (m_totEvtSize/4)%CIRC_BUFFER_SIZE;  
// Don't let the totEvtSize overflow; this is the integrated size over
// all events
            if(m_totEvtSize > 4*CIRC_BUFFER_SIZE) m_totEvtSize -= 4*CIRC_BUFFER_SIZE; 



// Make sure truncation bits are not set
            unsigned int Descriptor = *m_evtDescriptor; 
//            printf("Event Desc 0x%8.8x ...setting for no truncation\n",Descriptor);
            *m_evtDescriptor = (Descriptor & 0x3801ffff)| ((toWrite-EVT_HEADER_SIZE*4)/4) << 17;       
//            printf("New Event Descriptor 0x%8.8x....no truncation\n",*m_evtDescriptor);            
       }

// now write out the packet.
//       printf("Writing packetSize %i 0x%8.8x First word 0x%8.8x\n",packetSize,packetSize,*m_evtBuffer);
       if(writeEbf) swrite(m_fp, m_evtBuffer, packetSize/sizeof(*m_evtBuffer));

// If we are writing into the Tds save in a buffer, we must
// undo the swap if we wrote out the Ebf file
       if(TdsBuffer != 0) {       
         if(writeEbf) swap(m_evtBuffer,packetSize/sizeof(*m_evtBuffer));
         for (int i=0; i<packetSize/sizeof(*m_evtBuffer); i++) {
           TdsBuffer[i+offset] = m_evtBuffer[i];
           length=length+sizeof(*m_evtBuffer);
         }
         if(writeEbf) swap(m_evtBuffer,packetSize/sizeof(*m_evtBuffer));
         offset += packetSize/sizeof(*m_evtBuffer);
       }
             
    } while(truncate & !LdfFormat);


// OK everything is dumped to file clear buffer length.
    m_curEvtSize = 0;


   return (char *)m_evtBuffer;
}





    

/**
 *
 * @fn     int EbfOutput::close ()
 * @brief  Closes the file associated with this output stream.
 * @return Status
 *
 */
unsigned int EbfOutput::close ()
{
   int status;

   /* If already closed */
   if (m_fp == 0) return 0;
   
   if ((status = fclose (m_fp))) return status;

   //printf ("EbfOutput::close:      Have closed file status = %8.8x\n", status);
   m_fp = 0;
   
   return status;
}



/**
 *  
 *  @fn    EbfOutput::~EbfOutput ()
 *  @brief Closes the file associated with this stream (if it hasn't 
 *         already been closed) and frees the internal memory it 
 *         has allocated.
 */
EbfOutput::~EbfOutput ()
{
 
   /*
    | Must close the file if it is still opened and free the memory
    | associated with the event buffers, numbers and lengths. All
    | three of these blocks were allocated contigiously so just need
    | to free the memory associated with the top of this block. This
    | is the event buffer.
   */
  
   if (m_fp) close ();

   free (m_evtBuffer);


}





/**
 *
 *  @fn    void EbfOutput::print () const
 *  @brief Prints the event(s) currently in the output buffer to 
 *         standard output stream. 
 *
 */
void EbfOutput::print () const
{
    const unsigned int *evt = m_evtBuffer;
    int              totLen = m_curEvtSize;

    printf("************* Starting Dump of Ebf ***********\n");
 
    while (totLen > 0)
    {
        int elen = printOneEvt (evt);
        if (elen < (int)(18 * sizeof (struct LcbHeader)))
        {
            printf ("EbfOutput:print error, event too small %d\n", elen);
            return;
        }
        
        totLen -= elen;
    }
    
    return;
}



static int printOneEvt (const unsigned int *evt)
{
    

// Dump the event header words 
    for(int i=0; i<EVT_HEADER_SIZE-1; i++){
       printf("Event Header %i 0x%8.8x\n",i,evt[0]);   
       evt = (const unsigned int *)_advance(evt,4);
    }   

// Grab the event size (last of header words)
    int    evtLen = *evt++;
    int    remaining = evtLen - sizeof(int);
    printf("Event Length 0x%8.8x (%i) \n",evtLen,evtLen);
    int more = 1;

//    while (remaining > 0)
    while (more)    
    {
// This is the number of lcbWords (128 bits = 16 bytes)    
        int                clen = ((struct LcbHeader *)evt)->hdr_len & 0xffff;

// Convert lcbWords to bytes
        clen *= 16;        
        const unsigned int *nxt = _advance (evt, clen);


        if (clen < (int)sizeof (LcbHeader))
        {
            printf ("EbfOutput::printOneEvt error, contributor too small %d\n",
                    clen);
            return 0;
        }
        
        printOneContributor (evt, nxt, 71);
        remaining -= clen;
        more = (((struct LcbHeader *)evt)->hdr_len >> 31) & 0x1;
        evt        = nxt;
    }
    
// Rest of the data is padding;
    const unsigned int *nxt = _advance(evt, remaining);
    printf("Event Padding:\n");
    printData(evt,nxt,71);
    
    return evtLen;
}




static void printOneContributor (const unsigned int *beg, 
                                 const unsigned int *end,
                                 int         rightMargin)
/*
   DESCRIPTION
   -----------
   Prints the header and data from 1 contributor.

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
    static const char CNames[32][6] =
    {  
       "TWR 0",
       "TWR 1",
       "TWR 2",
       "TWR 3",
       "TWR 4",
       "TWR 5",
       "TWR 6",
       "TWR 7",
       "TWR 8",
       "TWR 9",
       "TWR A",
       "TWR B",
       "TWR C",
       "TWR D",
       "TWR E",
       "TWR F",
       "GEM  ",
       "ACD  ",
       "C x12",
       "C x13",
       "C x14",
       "C x15",
       "C x16",
       "C x17",
       "C x18",
       "C x19",
       "C x1A",
       "C x1B",
       "C x1C",
       "C x1D",
       "C x1E",
       "C x1F",
   };

    const struct LcbHeader *lcb = (const struct LcbHeader *)beg;
    int                     cid = (lcb->hdr_len >> 17) & 0x1f;

    /* Print the LCB Header */
    printf ("\n%s Hdr:%8.8x %8.8x\n", CNames[cid], beg[0], beg[1]);
 
    /* Print the Data block */
    printData ((const unsigned int *)(lcb + 1), end, rightMargin);
    

    
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
