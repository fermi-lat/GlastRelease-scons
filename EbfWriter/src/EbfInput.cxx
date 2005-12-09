#include <stdio.h>
#include <stdlib.h>

#include "EbfAcdData.h"
#include "EbfCalData.h"
#include "EbfTkrData.h"
#include "EbfGemData.h"
#include "EbfInput.h"
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
static inline unsigned int      sread (FILE         *fp,
                                        unsigned int *wrds,
                                        int          nwrds);

//static inline void              parse(unsigned int *evt);


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
       while (--npad >= 0) *dst++ = 0x55;
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
       while (--npad >= 0) *dst++ = 0xFF;
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

   /* Compute the parity over the 6 bits of the source and destination */
   parity       =  sourceId ^ destinationId;
   parity       =  parity ^ (parity >> 3);
   parity       = (parity ^ (parity >> 1) ^ (parity >> 2)) & 1;
   /* The bit at 15 is the response bit which indicates that additional
   / contributions follow.  We set this for all an then zero the bit in
   / the last contribution when we know what it is.
   */
   return (1 << 15) | (destinationId << 9) | (LATp & 0x3) << 7 | (sourceId << 1) | (parity ^ 1);
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
  int trigAck   = 0;
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

   hdr->hdr_len = (headerWord (sourceId, destinationId) << 16) | lcbWords;
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
//     printf("Before Swap 0x%8.8x \n",tmp);
     tmp = ((tmp & 0x000000ff) << 24) |
       ((tmp & 0x0000ff00) <<  8) |
       ((tmp & 0x00ff0000) >>  8) |
       ((tmp & 0xff000000) >> 24); 
//     printf("After Swap 0x%8.8x \n",tmp);  
     *wrds++ = tmp;
   }
   
   return;
}






 unsigned int sread  (FILE           *fp,
                            unsigned int *wrds,
                            int          nwrds)
/*
    DESCRIPTION
    -----------
    Reads the specified number of 32-bit words to the specified output
    stream. If necessary, the 32-bit words are byte-swapped into 
    big-endian format. Note, if the buffer is byte-swapped, it is done
    in place.

    PARAMETERS
    ----------
            fp: The file handle of the output stream.

          wrds: The 32-bit words to Read in.

         nwrds: The number of 32-bit words to Read in.

   RETURNS
   -------
   The return value of .
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
    unsigned int val = fread (wrds, sizeof (*wrds), nwrds, fp);
    swap (wrds, nwrds);
    return val;
}






/**
 *
 *  @fn      EbfInput::open (const char    *fileName,
 *                            unsigned int maxEvtSize)
 *  @brief Initializes the EbfInput be opening the named file and allocating
 *         enough memory to hold at least 1 maximally sized event.
 *
 *  @param fileName   The name of the file to open.
 *  @param maxEvtSize The maximum size of a given event.
 *
 */
int EbfInput::open (const char     *fileName,
                     unsigned int  maxEvtSize,
                     bool ldfFormat)
{
   unsigned char           *ptr;

   printf("EbfInput::open file %s\n",fileName);

   /* Initialize all the data members to something beign */   
   m_fp         = 0;
   m_maxEvtSize = maxEvtSize;
   m_curEvtSize = 0;
   m_totEvtSize = 0;
   m_numEvtsIn  = 0;
   m_numEvtsOut = 0;
   m_evtBuffer  = 0;
   m_evtHeader  = m_evtBuffer;
   m_nReadEvt   = 0;

   
   /* Allocate a buffer big enough to hold a maximally sized event */
   ptr = (unsigned char *)malloc (maxEvtSize);
   if (ptr == 0) { return -1;  }

   
   /* Open the output file */
   char file[120];
   ldfFormat ? sprintf(file,"%s.ldf",fileName) : sprintf(file,"%s.ebf",fileName);
   m_fp = fopen (file, "rb");
   if (!m_fp)
   {
       /* Error in opening the output file.. */
     free (ptr);
     printf("EbfInput: Error openning the input file\n");
     return -2;
   }
   

   m_evtBuffer = (unsigned int *)ptr;
   m_evtHeader = m_evtBuffer;


   /* Report the name of the output file, maximum event size */
   /*
     printf ("EbfInput::initialize: FileName: %s\n"
     "                 Max Event Size: %d\n",
     fileName,
     maxEvtSize);
   */
   return 0;
}


/**
 *
 *  @fn     unsigned int EbfInput::read ()
 *  @brief  Writes the current event to the output stream.
 *  @return The number of 32 bit words written
 *
 */
unsigned int EbfInput::read (bool ldfFormat)
{
    
   int toRead = 0;
   bool debug = true;
   
   if(!ldfFormat) {
       if(debug) printf("EbfInput: Reading from File...EBF format\n");
    // Start by Reading in the Event Header
       unsigned int sizeRead = sread(m_fp, m_evtHeader, EVT_HEADER_SIZE);
       m_evtBuffer = m_evtHeader;
       for(int i=0; i<EVT_HEADER_SIZE; i++) {
          if(debug) printf("Word %i = 0x%8.8x \n",i,m_evtBuffer[0]);
          if(i==EVT_HEADER_SIZE-1) toRead = (m_evtBuffer[0]>>17)&0x7ff;
          m_evtBuffer = (unsigned int *)_advance(m_evtBuffer,4);
       }   

       if(debug) printf("Event Length to Read 0x%8.8x \n",toRead);
       if(debug) sizeRead = sread (m_fp, m_evtBuffer, toRead);
   } else {

// Format is LDF form
       if(debug) printf("EbfInput: Reading from File...LDF format: Event %i\n",m_nReadEvt);
       m_nReadEvt++;
       bool ebfContrib = false;
       unsigned int contribType = 0;
       
       while(!ebfContrib) {
// Start by Reading in the Event Header
         unsigned int sizeRead = sread(m_fp, m_evtHeader, LDF_HEADER_SIZE);
         if(debug) printf("Size Read %i  Event %i \n",sizeRead,m_nReadEvt);

//         m_evtHeader = m_evtBuffer;
         m_evtBuffer = m_evtHeader;       
         for(int i=0; i<LDF_HEADER_SIZE; i++) {
            if(debug) printf("Word %i = 0x%8.8x \n",i,m_evtBuffer[0]);
            if(i==LDF_HEADER_SIZE-1) toRead = m_evtBuffer[0]-8;
            if(i==LDF_HEADER_SIZE-2) contribType = m_evtBuffer[0] & 0xfffff; 
            m_evtBuffer = (unsigned int *)_advance(m_evtBuffer,4);
         }   

         if(debug) printf("Event Length to Read 0x%8.8x ContribType 0x%8.8x\n",toRead,contribType);
         sizeRead = sread (m_fp, m_evtBuffer, toRead/4);
       
         if(contribType == 0xF0010) {
            if(debug) printf("EbfInput:Read  EBF Contributor Found\n");
            ebfContrib = true;
//            m_nReadEvt++;
         } else {
            if(debug) printf("EbfInput:Read  LDF Contributor not EBF:  skip\n");
         }
       } //while looking for EBF contributor 
   
   }    
   
   return toRead;
}


/**
 *
 *  @fn     unsigned int EbfInput::parse ()
 *  @brief  Writes the current event to the output stream.
 *  @return The number of 32 bit words written
 *
 */ 
void EbfInput::parse(   EbfTkrData *tkr, 
                        EbfCalData *cal,
                        EbfGemData *gem, 
                        EbfAcdData *acd,
                        EbfCalConstants *calCon)
{


    unsigned int *evt = m_evtBuffer++;
    bool more = true;
    bool debug = true;

// Cycle through the Contributions
   while (more) { 

// Read Contribution Header
      int contHeader = *evt;
      if(debug) printf("Contribution Header 0x%8.8x\n",contHeader);

// Get the number of lcb Words      
      int lcbWords = contHeader & 0xffff;
      int sourceID = (contHeader >> 17) & 0x3f;
      ((contHeader&0x80000000) >> 31)==1 ? more=true : more=false;
      if(debug) printf("SourceID %i #lcbWords %i More %i\n",sourceID,lcbWords,more);
      
// Sort Out the contribution
      if(sourceID>=EbfContributorId::TWR  && sourceID<EbfContributorId::TWR+16) {

// TEM Contribution
         int itower = sourceID-EbfContributorId::TWR;
         
         if(debug) printf("Found TEM Contribution for Tower %i:\n",itower);
         cal->parseInput(evt,itower,lcbWords,calCon);
         tkr->parseInput(evt,itower,lcbWords);

         evt = (unsigned int *)_advance(evt,lcbWords*16);

      }else if(sourceID==EbfContributorId::GEM) {

// GEM Contribution
         if(debug) printf("Found GEM Contribution: \n");
         gem->parseInput(evt,lcbWords);      
         evt = (unsigned int *)_advance(evt,lcbWords*16);
         
      }else if(sourceID==EbfContributorId::ACD) {

// ACD Contribution
         if(debug) printf("Found ACD Contribution:\n");
         acd->parseInput(evt,lcbWords);
         evt = (unsigned int *)_advance(evt,lcbWords*16);
      } else {

// Some other contribution advance the event pointer
         if(debug) printf("Found nonData Contribution:\n");
         evt = (unsigned int *)_advance(evt,lcbWords*16);      
      }          
      
   }
      
   return;
}



/**
 *
 * @fn     int EbfInput::close ()
 * @brief  Closes the file associated with this output stream.
 * @return Status
 *
 */
unsigned int EbfInput::close ()
{

   int status = 0;

//   printf("Closing Input EbfFiles\n");
   // If already closed 
   if (m_fp == 0) return 0;
   
   if ((status = fclose (m_fp))) return status;

   //printf ("EbfInput::close:      Have closed file status = %8.8x\n", status);
   m_fp = 0;
   
   return status;
}



/**
 *  
 *  @fn    EbfInput::~EbfInput ()
 *  @brief Closes the file associated with this stream (if it hasn't 
 *         already been closed) and frees the internal memory it 
 *         has allocated.
 */
EbfInput::~EbfInput ()
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
//   printf("In the Destructor for EbfInput\n");
   if (m_fp) close ();

//   free (m_evtBuffer);


}





/**
 *
 *  @fn    void EbfInput::print () const
 *  @brief Prints the event(s) currently in the output buffer to 
 *         standard output stream. 
 *
 */
void EbfInput::print () const
{
    const unsigned int *evt = m_evtBuffer;
    int              totLen = m_curEvtSize;

    printf("************* Starting Dump of Ebf ***********\n");
 
    while (totLen > 0)
    {
        int elen = printOneEvt (evt);
        if (elen < (int)(18 * sizeof (struct LcbHeader)))
        {
            printf ("EbfInput:print error, event too small %d\n", elen);
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
            printf ("EbfInput::printOneEvt error, contributor too small %d\n",
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
