#include <stdio.h>
#include <stdlib.h>

#include "EbfAcdData.h"
#include "EbfCalData.h"
#include "EbfTkrData.h"
#include "EbfGltData.h"
#include "EbfOutput.h"

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


/* Define a dummy EPU to act as the destination */
#define EPU 0x18


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
static inline unsigned int summaryWord (unsigned int sequence);
static inline unsigned int         pad (const unsigned int *beg,
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
       while (--npad >= 0) *dst++ = 0x55;
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

   sourceId      &= 0x1f;
   destinationId &= 0x1f;

   /* Compute the parity over the 6 bits of the source and destination */
   parity       =  sourceId ^ destinationId;
   parity       =  parity ^ (parity >> 3);
   parity       = (parity ^ (parity >> 1) ^ (parity >> 2)) & 1;
   return (destinationId << 9) | (sourceId << 1) | (parity ^ 1);
}




static unsigned int summaryWord (unsigned int sequence)
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
  return  ((sequence & 0x3) << 29) | (((sequence >> 2) & 0x3fff) << 1);
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

   hdr->hdr_len = (headerWord (sourceId, destinationId) << 16) | length;
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
    return fwrite (wrds, sizeof (*wrds), nwrds, fp);
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


   /* Initialize all the data members to something beign */   
   m_fp         = 0;
   m_maxEvtSize = maxEvtSize;
   m_curEvtSize = 0;
   m_totEvtSize = 0;
   m_numEvtsIn  = 0;
   m_numEvtsOut = 0;
   m_evtBuffer  = 0;

   
   /* Allocate a buffer big enough to hold a maximally sized event */
   ptr = (unsigned char *)malloc (maxEvtSize);
   if (ptr == 0) { return -1;  }

   
   /* Open the output file */
   //m_fp = fopen (fileName, "wb");
   //if (!m_fp)
   //{
       /* Error in opening the output file.. */
   //free (ptr);
   //return -2;
   //}
   

   m_evtBuffer = (unsigned int *)ptr;


   /* Report the name of the output file, maximum event size */
   /*
     printf ("EbfOutput::initialize: FileName: %s\n"
     "                 Max Event Size: %d\n",
     fileName,
     maxEvtSize);
   */
   return 0;
}





/**
 *
 *  @fn   unsigned int EbfOutput::format (const EbfAcdData *acd,
 *                                        const EbfCalData *cal,
 *                                        const EbfTkrData *tkr,
 *                                        const EbfGltData *glt)
 *  @brief Converts the detector data blocks into EBF format.
 *  @param acd  The ACD data block.
 *  @param cal  The CAL data block, contains all 16 towers.
 *  @param tkr  The TKR data block, contains all 16 towers.
 *  @param glt  The GLT data block.
 *  @return     The number of 32 bit words needed to format the event
 *
 */
unsigned int EbfOutput::format (const EbfAcdData *acd,
                                const EbfCalData *cal,
                                const EbfTkrData *tkr,
                                const EbfGltData *glt)
{
    int print = m_print;
    

   /*  
    | This convention matches the event display's numbering convention,
    | ie. the first event out is labeled #1. This number will be tacked
    | on to the end of the GLT data (hidden away).
    */
   unsigned int numEvtsIn = m_numEvtsIn + 1;
   
   
   /* Count the number of events seen */
   m_numEvtsIn = numEvtsIn;
 

   /* Check if this event is triggered */
   if (!glt->isTriggered ())
   {
       //printf (" WARNING: Event in = %8d not triggered\n", numEvtsIn);
       return 0;
   }

   
   unsigned int numEvtsOut = m_numEvtsOut;


   /* Count the number of events output */   
   m_numEvtsOut = numEvtsOut + 1;


   /* Event is triggered, format it to the output buffer */
   unsigned int  summary = summaryWord (numEvtsOut); // Event sequence #
   unsigned int  *evtBeg = m_evtBuffer;              // Event begin
   unsigned int     *beg = evtBeg + 1;               // Contributor begin
   unsigned int     *dst;
   

   /* Pack up the GLT contribution */
   dst    = _advance (beg, sizeof (LcbHeader));
   dst    = glt->format (dst);
   dst    = complete (beg, dst, EbfContributorId::GLT, EPU, summary);
   /*
     if (print)
     {
     printf ("GLT data:\n"); 
     printData (beg, dst, 72);
     }
   */
   


   /* Fill in the Tower Contributions */
   for (unsigned int itower = 0; itower < EbfCalData::NumTowers; itower++)
   {
       beg = dst;
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
   

   /* Fill in the event length and event number of this event */
   unsigned int length = (unsigned char *)dst - (unsigned char *)evtBeg;
  *evtBeg        = length;
   m_curEvtSize  = length;
   m_totEvtSize += length;

      
   /* Return the length */
   return length;
}





/**
 *
 *  @fn     unsigned int EbfOutput::write ()
 *  @brief  Writes the current event to the output stream.
 *  @return The number of 32 bit words written
 *
 */
char * EbfOutput::write (unsigned int &length)
{
   size_t toWrite = m_curEvtSize;

   m_curEvtSize = 0;
   swap(m_evtBuffer,toWrite/sizeof(*m_evtBuffer));
   length=toWrite;
   //swrite (m_fp, m_evtBuffer, toWrite/sizeof(*m_evtBuffer));
   //return swrite (m_fp, m_evtBuffer, toWrite/sizeof(*m_evtBuffer));
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
   /* If have anything left in the event buffer, write it out */
  //if (m_curEvtSize) write ();
 
 
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
    int    evtLen = *evt++;
    int remaining = evtLen - sizeof (int);
    
    while (remaining > 0)
    {
        int                clen = ((struct LcbHeader *)evt)->hdr_len & 0xffff;
        const unsigned int *nxt = _advance (evt, clen);

        if (clen < (int)sizeof (LcbHeader))
        {
            printf ("EbfOutput::printOneEvt error, contributor too small %d\n",
                    clen);
            return 0;
        }
        
        printOneContributor (evt, nxt, 72);
        remaining -= clen;
        evt        = nxt;
    }
    
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
       "GLT  ",
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
