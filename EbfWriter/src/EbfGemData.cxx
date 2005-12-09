#include <stdio.h>

#include "EbfGemData.h"
#include "EbfAcdData.h"
#include "EbfCalData.h"
#include "EbfTkrData.h"
#include "EbfGemCounters.h"


/*
 | This set of routines uses data from the ACD, TKR and CAL to form
 | the GEM's contribution to the event data.
 |
 | There are 2 non-trivial calculation in this file. The first is the
 | calculation of the ACD throttle tower mask and the second  is the
 | handling of time.
 |
 | ACD THROTTLE MASK CALCULATION
 | One may ask why is the ACD throttle formed here rather than in
 | the ACD code, as in the other systems. The reason is that the
 | logic forming the ACD throttle is on the GEM not the ACD. The
 | tracker 3-in-a-row and the CAL HI/LO are formed on the TKR/CAL
 | TEM board, so the code doing that is there. The ENABLES/DISABLES
 | for which CAL HI/LO signals go into the formation are also on the
 | individual TEMS. The ENABLES/DISABLES for which trigger signals
 | are used are in the GEM. The CNO trigger is really a collection
 | of 12 signals, one from each board. It is the responsibility of
 | the ACD code to generate the 12 signals. It is the responsibilty
 | of the GEM code to turn them into a trigger.
 |
 |
 | !!! IMPROVEMENT OPPORTUNITY !!!
 | -------------------------------
 |
 | TIME HANDLING and GEM Modeling
 | These registers and their underlying supporting counters are part
 | of the GEM. However, one could make a case that the should not
 | be modelled here, only the information they contain should be 
 | extracted here. This is much more in keeping with the philosophy
 | of the other subsystems, ie the ACD, CAL and TKR. However, GLEAM,
 | in the state that I received it did not model the GEM (it only
 | had some sub-portion of the trigger decision logic implemented)
 | in a separate class. I did not feel I was in a position to do this
 | correctly, but a future major improvement to GLEAM should include
 | a more complete model of the GEM implemented in an architecturally
 | consistent way.
 |
*/ 
static unsigned int formThrottle       (unsigned int                  xz,
                                        unsigned int                  yz,
                                        unsigned int                  xy);
static unsigned int formThrottleByFace (int                        tiles, 
                                        const unsigned short int *towers);



/**
 *
 *  @fn     void fill (unsigned int       mc_number,
 *                     double               mc_time,
 *                     EbfGemCounters *lat_counters,
 *                     const EbfAcdData        *acd,
 *                     const EbfTkrData        *tkr,
 *                     const EbfCalData        *cal)
 *  @brief  Fills  the GEM Record and serves as the basis of the 
 *          trigger decision.
 *
 *  @param  mc_number  The identifying GLEAM event number. This
 *                     is squirreled away in the GEM so that it
 *                     can be put into a hidden location of the
 *                     GEM output record.
 *  @param  mc_time    The GLEAM event generation time. This is
 *                     used to simulate the LAT time registers
 *                     counters.
 *  @param lat_counter Pointer to a class that models the counters
 *                     and their supporting registers on the GEM.
 *  @param acd         Pointer to a class storing the ACD data in
 *                     a format that can be easily used to form
 *                     the GEM information.
 *  @param tkr         Pointer to a class storing the TKR data in
 *                     a format that can be easily used to form
 *                     the GEM information.
 *  @param cal         Pointer to a class storing the CAL data in
 *                     a format that can be easily used to form
 *                     the GEM information.
 */
void EbfGemData::fill (Event::EventHeader *header,
                       unsigned int       mc_number,
                       double               mc_time,
                       EbfGemCounters *lat_counters,
                       const EbfAcdData        *acd,
                       const EbfTkrData        *tkr,
                       const EbfCalData        *cal)
{
    unsigned int            xz;
    unsigned int            yz;
    unsigned int            xy;
    unsigned int            ru;
    unsigned short threeInARow;
    unsigned short   throttled;
    unsigned short       calHi;
    unsigned short       calLo;
    unsigned short         cno;
    unsigned short      reqVec;  


    xz = acd->vetoesXZ ();
    yz = acd->vetoesYZ ();
    xy = acd->vetoesXY ();
    ru = acd->vetoesRU ();
    
/*    printf("xz: 0x%8.8x\n",xz);
    printf("yz: 0x%8.8x\n",yz);
    printf("xy: 0x%8.8x\n",xy);
    printf("ru: 0x%8.8x\n",ru);
*/
    /*
     | Form the following quantities:
     |    0. Throttled towers by the ACD            (per tower) 
     |    1. Track 3-in-a-row                       (per tower)
     |    2. CAL HI                                 (per tower)
     |    3. CAL LO                                 (per tower)
     |    4. ACD CNO                                (global   )
     |
    */
    threeInARow = tkr->threeInARow ();
    throttled   = formThrottle (xz, yz, xy);
    calHi       = cal->hiTrigger ();
    calLo       = cal->loTrigger ();
    cno         = acd->cno ();

/*
   printf("Ebf Alg Threshold \n");
   printf("Hit Towers    0x%8.8x\n",threeInARow);
   printf("ACD Mask xy   0x%8.8x\n",xy);
   printf("ACD Mask yz   0x%8.8x\n",yz);
   printf("ACD Mask xz   0x%8.8x\n",xz);
   printf("--------------------------\n");
*/

    /* Collapse the triggering information into the trigger request vector */
    reqVec      = ((threeInARow & throttled)  ? (1 << TKR_NOT_VETOED) : 0)
                | ((threeInARow)              ? (1 << TKR)            : 0)
                | ((calHi)                    ? (1 << CAL_HI)         : 0)
                | ((calLo)                    ? (1 << CAL_LO)         : 0)
                | ((cno)                      ? (1 << CNO)            : 0);


    unsigned int dif = (reqVec&0xff^( (header->trigger()&0xff00)>>8) ); 
//    if(dif>0) printf("Cond. Sum Difference:  EBF 0x%2.2x    Header 0x%2.2x  HeaderGlt 0x%2.2x  XOR 0x%2.2x\n",reqVec,header->trigger(),(header->trigger()&0xff00)>>8,dif);

    /* Copy the trigger information to the GEM */
    m_thrTkr    = ((throttled & 0xffff)<< 16) | (threeInARow & 0xffff) ;
    m_calHiLo   = ((calHi & 0xffff)    << 16) | (calLo & 0xffff);
    m_cnoReqvec = ((reqVec & 0xff)     << 16) | (cno & 0xfff);


    /* Copy the ACD hit data from the ACD record to the GEM */
    m_tiles[XZ] = xz;
    m_tiles[YZ] = yz;
    m_tiles[XY] = xy;
    m_tiles[RU] = ru;

    /* Keep track of the Monte Carlo number and the event time */
    m_number    = mc_number;
    m_time      = mc_time;


    /* 
     | Fill in the GEM busy, prescale and sent counters.
     |  
     | !!! KLUDGE !!!
     | --------------
     | These quantities need to be emulated better.
     |
     | SENT COUNTER
     | Independent of whether this event generated a trigger, it
     | did cause a window turn. Unfortunately there are lot of
     | other window turns that GLEAM has no way of modeling. There
     | are many other 'events' that can generate a window turn, but
     | not result in a trigger.
     |
     | BUSY AND PRESCALE COUNTERS
     | The busy and prescaled are not being modeled at all. The 
     | window turn counter (sent) has, at least, a very crude model
     | of incrementing on every generated event.
    */
    m_busy       = lat_counters->busy () & 0xffffff;
    m_prescaled  = lat_counters->prescale ()  & 0xffffff;
    m_sent       = lat_counters->increment_sent ()  & 0xffffff;


    /*
     | TIME REGISTER AND COUNTERS
     | This section fills in the GEM's time registers that eventually
     | go into the GEM's contribution to the event data. These are
     |
     |    a) The time of event in terms of LAT system clock ticks
     |    b) The 1 PPS register
     |    c) The live time
     |
     | The first two, the event time and 1 PPS register are sourced
     | from the same underlying counter, the 25 bit LAT system clock
     | counter. The event time is simply the value of this counter
     | at the time the trigger is declared. The 1 PPS time is value
     | of this counter at the time of the 1 PPS hack. The 1 PPS 
     | register also contains, in the upper 7 bits, an value that
     | increments with every 1 PPS time hack. This is used to 
     | correlate the 1 PPS time in this register with a 1 PPS message.
     |
     | The live time register is simply the LAT system clock gated
     | by the DAQ LIVE signal (or if you perfer NOT DAQ BUSY.
     |
    */


    /* Convert the total elapsed time in system clock ticks */
    unsigned int lat_ticks     = lat_counters->ticks (mc_time);


    /* This gives the value of the PPS counter at the last whole second */
    unsigned int whole_seconds = (unsigned int)mc_time;
    unsigned int pps_ticks     = lat_counters->ticks (whole_seconds);
       
 
    /*
     | Pack the event time and 1 PPS registers.The event time is low
     | 25 bits of the lat time in ticks. Similarly the 1 PPS register
     | contains the low 25 bits of the lat time in ticks at the last
     | whole second. Remember the upper 7 bits contain in number that
     | increments with every 1 PPS time hack. Technically this does
     | have to be correlated with 1 second passing (for example, if
     | the 1 PPS time hack is missed), but for the purposes here this
     | number will be treated as the lower 7 bits of the number of
     | whole elapsed seconds. Since the number of elapsed seconds
     | may exceed 128 seconds, the reading software is processing this
     | information is expected to notice that the number of seconds 
     | goes from 127 => 0, and interpret this as a signal to add 128
     | seconds to the total lat time in ticks. 
    */
    m_evttime  = lat_ticks & 0x1ffffff;
    m_ppstime  = pps_ticks & 0x1ffffff;
    m_ppstime |= (whole_seconds & 0x7f) << 25;
    m_deltaTime = m_evttime - m_lastEvtTime;
    m_lastEvtTime = m_evttime;    

//    printf("Event Time %f gemTime 0x%8.8x ppstime 0x%8.8x seconds 0x%8.8x \n",
//            mc_time,m_evttime,m_ppstime,whole_seconds); 

    /*
     | LIVE TIME COUNTER 
     | This computes the dead time contribution of this event. The model
     | used is extremely crude, each event, independent of size or
     | anyother conditions contributions 20 useconds of dead time. 
     |
     | !!! KLUDGE !!!
     | --------------
     | A more sophisticated model based on component contributions
     | and past history. The LAT DAQ is heavily buffered so paths may
     | be busy from previous activities is needed. Stated more correctly,
     | this event may cause these paths to become busy thus blocking
     | future triggers for some time in the future. This is really hard
     | to do correctly. One needs to take into account all the state of
     | all buffering FIFOs, (TKR, TEM, EB) and even the amount of time
     | it takes an EPU to process an event. Likely the best one can do
     | is come up with an empirical formula that roughly mimics this.
     | Good luck!!!
    */
    unsigned int dead_ticks = lat_counters->declare_dead (.000020);
    
    /* Fill in the live time counter */
    m_livetime = (lat_ticks - dead_ticks) & 0xffffff;
    
    return;
}




/**
 *
 *  @fn     unsigned int *EbfGemData::format (unsigned int *dst)
 *  @brief  Converts the GEM data into Event Builder Format.
 *
 *  @param  dst The destination array.
 *  @return     Pointer to the next available word in the output array.
 * 
 *  It is the responsibility of the caller to ensure that enough room
 *  for a maximum record exists in the output array.  
 *
 */
unsigned int *EbfGemData::format (unsigned int *dst) const
{
   dst[0]  = m_thrTkr;
   dst[1]  = m_calHiLo;
   dst[2]  = m_cnoReqvec;
   dst[3]  = m_tiles[XZ];
   dst[4]  = m_tiles[YZ];
   dst[5]  = m_tiles[XY];
   dst[6]  = m_tiles[RU];
   dst[7]  = m_livetime;
   dst[8]  = m_prescaled;
   dst[9]  = m_busy;   
   dst[10] = m_sent;
   dst[11] = m_evttime;
   dst[12] = m_ppstime;
   dst[13] = m_deltaTime;

   
   return &dst[14];
}


void EbfGemData::parseInput(unsigned int *contrib, unsigned int lcbWords) const
{

   unsigned int *data = contrib;
   bool debug = false;
   if(debug) {
      printf("GEM Contribution:\n");
      for(int i=0; i<lcbWords; i++) {
         for(int j=0; j<4; j++) {
             printf(" 0x%8.8x ",*data++);
         }
         printf("\n");
      }
      printf("\n");
   // reset pointer
      data = contrib;
   }
   
   return;
}



/**
 *
 *  @fn     void print ()
 *  @brief  Diagnostic print routine 
 * 
 *  Routine to print the GEM data for this event. This is mainly used
 *  used to debug the fill and format routines,
 *
 */
void EbfGemData::print() const
{
   std::printf ("GEM    ThrTkr = %8.8x\n"
                "      CalHiLo = %8.8x\n"
                "    CnoReqvec = %8.8x\n"
                "     XZ Tiles = %8.8x     YZ Tiles = %8.8x\n"
                "     XY Tiles = %8.8x     RU Tiles = %8.8x\n"
                "     Livetime = %8.8x   Busy Count = %8.8x\n"
                "    Prescaled = %8.8x   Sent Count = %8.8x\n"
                "      EvtTime = %8.8x      PpsTime = %8.8x\n",
                m_thrTkr,
                m_calHiLo,
                m_cnoReqvec,
                m_tiles[XZ],
                m_tiles[YZ],
                m_tiles[XY],
                m_tiles[RU],
                m_livetime,
                m_busy,
                m_prescaled,
                m_sent,
                m_evttime,
                m_ppstime);
       
   return;
}





/*
 |  ACD THROTTLE COMPUTATION
 |  ------------------------
 |  The following internal routines are used to compute which towers are
 |  shadowed by the collection of ACD hits in this event. Note that this
 |  calculation does not reference the TKR data at all. It is strictly
 |  a mapping of a set of struck ACD tiles to a set of towers. 
 |
 |
 |  !!! IMPROVEMENT OPPORTUNITY!!!
 |  ------------------------------
 |  A data statement in the code determines the mapping of a struck ACD
 |  tile to a set of towers. A better implementation would promote this
 |  information up to the user interface level. A intermediate solution
 |  would be provide a palette of a couple of realistic mappings that
 |  is user selectable.
 |
 */



/**
 *
 *  @def    TWR(_n) 
 *  @brief  Transforms a tower number, specified as 0 - f, into its
 *          corresponding bit list. Tower 0 = LSB.
 *
 *  @param  _n The tower number, specified as 0 - f. Note that no '0x'
 *             is used to denote a hex number.
 *  @return     An integer with the bit corresponding to the tower set.
 *              The convention is LSB = Tower 0.
 *
**/
#define TWR(_n)   (1 << 0x ## _n)


/**
 *  @def    IGNORE
 *  @brief  Defines the list of ACD side tiles to ignore when forming
 *          the throttle
 *
 *   Currently, the lower two rows (the big tile and the 5 tiles in the
 *   immediately above it) are ignored.
 *
**/
#define IGNORE    ((1<<10) | (1<<11) | (1<<12) | (1<<13) | (1<<14) | (1<<15))


/**
 *
 *  @def    XZ_IGNORE
 *  @brief  Defines the ignore list for both XZ+ and XZ- plane
 * 
**/
#define XZ_IGNORE ((IGNORE<<16) | IGNORE)


/**
 *
 *  @def    YZ_IGNORE
 *  @brief  Defines the ignore list for both XZ+ and XZ- plane
 * 
**/
#define YZ_IGNORE ((IGNORE<<16) | IGNORE)


static unsigned int formThrottle (unsigned int xz,
                                  unsigned int yz,
                                  unsigned int xy)
/*
   DESCRIPTION
   -----------
   Computes which towers are considered shadowed according to the which
   tiles are struck. This is accomplished by examining each bit in the
   list of struck tiles. For each tile that is hit, a mask of towers
   considered shadowed by that tile is accumulated into a 16 bit mask,
   representing the shadowed towers.

   PARAMETERS
   ----------
           xz: The list of struck tiles in the ACD XZ planes 
               The upper 16 bits are associated with the XZ plane
               intersecting the +Y axis. 
               The lower 16 bits are associated with the XZ plane
               intersecting the -Y axis.

           yz: The list of struck tiles in the ACD YZ planes 
               The upper 16 bits are associated with the YZ plane
               intersecting the +X axis. 
               The lower 16 bits are associated with the YZ plane
               intersecting the -X axis.

           xy: The list of struck tiles in the ACD XY (the TOP) plane 

   RETURNS
   -------
   The list of shadowed towers. This list is right justified, with
   tower 0 represented by bit 0
*/         
{
   struct Shadowed
   {
       unsigned short int xz[32]; /* XZ face tiles shadowed towers */
       unsigned short int yz[32]; /* YZ face tiles shadowed towers */
       unsigned short int xy[25]; /* XY face tiles shadowed towers */
   };
      

   /*
    | This defines the mapping of an ACD tile to a set of shadowed towers.
    | Note that this table defines the shadowed towers even for tiles
    | that are potentially ignored. I wanted to keep the concepts of 
    | the set mapping and the ignore lists orthogonal.
    |
    | A word of warning. This structure is laid out as 3 arrays, one
    | for each of the XZ, YZ and XY faces. The ordering of the elements
    | within the array is in the order of MSB to LSB as the tiles appear
    | in the ACD struck tile bit map lists.
   */
   static const struct Shadowed ShadowedTowers =
   {
       /* --------------------------------------------- */
       /* Shadowed towers by XZ ACD tiles               */
       /* --------------------------------------------- */
       {
          /* ------------------------------------------ */
          /* XZ + Tiles                                 */
          /* ------------------------------------------ */
          TWR(C) | TWR(D) | TWR(E) | TWR(F), /* Tile 15 */

          TWR(F),                            /* Tile 14 */
          TWR(E) | TWR(F),                   /* Tile 13 */
          TWR(D) | TWR(E),                   /* Tile 12 */
          TWR(C) | TWR(D),                   /* Tile 11 */
          TWR(C),                            /* Tile 10 */
          
          TWR(F),                            /* Tile  9 */
          TWR(E) | TWR(F),                   /* Tile  8 */
          TWR(D) | TWR(E),                   /* Tile  7 */
          TWR(C) | TWR(D),                   /* Tile  6 */
          TWR(C),                            /* Tile  5 */
          
          TWR(F),                            /* Tile  4 */
          TWR(E) | TWR(F),                   /* Tile  3 */
          TWR(D) | TWR(E),                   /* Tile  2 */
          TWR(C) | TWR(D),                   /* Tile  1 */
          TWR(C),                            /* Tile  0 */
          

          /* ------------------------------------------ */
          /* XZ - Tiles                                 */
          /* ------------------------------------------ */
          TWR(0) | TWR(1) | TWR(2) | TWR(3), /* Tile 15 */
          
          TWR(3),                            /* Tile 14 */
          TWR(2) | TWR(3),                   /* Tile 13 */
          TWR(1) | TWR(2),                   /* Tile 12 */
          TWR(0) | TWR(1),                   /* Tile 11 */
          TWR(0),                            /* Tile 10 */
          
          TWR(3),                            /* Tile  9 */
          TWR(2) | TWR(3),                   /* Tile  8 */
          TWR(1) | TWR(2),                   /* Tile  7 */
          TWR(0) | TWR(1),                   /* Tile  6 */
          TWR(0),                            /* Tile  5 */
          
          TWR(3),                            /* Tile  4 */
          TWR(2) | TWR(3),                   /* Tile  3 */
          TWR(1) | TWR(2),                   /* Tile  2 */
          TWR(0) | TWR(1),                   /* Tile  1 */
          TWR(0)                             /* Tile  0 */
       },
       /* --------------------------------------------- */
       /* Shadowed towers by YZ ACD tiles               */
       /* --------------------------------------------- */       
       {
          /* ------------------------------------------ */
          /* YZ + tiles                                 */
          /* ------------------------------------------ */
          TWR(3) | TWR(7) | TWR(B) | TWR(C), /* Tile 15 */

          TWR(F),                            /* Tile 14 */
          TWR(B) | TWR(F),                   /* Tile 13 */
          TWR(7) | TWR(B),                   /* Tile 12 */
          TWR(3) | TWR(7),                   /* Tile 11 */
          TWR(3),                            /* Tile 10 */
      
          TWR(F),                            /* Tile  9 */
          TWR(B) | TWR(F),                   /* Tile  8 */
          TWR(7) | TWR(B),                   /* Tile  7 */
          TWR(3) | TWR(7),                   /* Tile  6 */
          TWR(3),                            /* Tile  5 */
          
          TWR(F),                            /* Tile  4 */
          TWR(B) | TWR(F),                   /* Tile  3 */
          TWR(7) | TWR(B),                   /* Tile  2 */
          TWR(3) | TWR(7),                   /* Tile  1 */
          TWR(3),                            /* Tile  0 */


          /* ------------------------------------------ */
          /* YZ - tiles                                 */
          /* ------------------------------------------ */
          TWR(0) | TWR(4) | TWR(8) | TWR(C), /* Tile 15 */
          
          TWR(C),                            /* Tile 14 */
          TWR(8) | TWR(C),                   /* Tile 13 */
          TWR(4) | TWR(8),                   /* Tile 12 */
          TWR(0) | TWR(4),                   /* Tile 11 */
          TWR(0),                            /* Tile 10 */
          
          TWR(C),                            /* Tile  9 */
          TWR(8) | TWR(C),                   /* Tile  8 */
          TWR(4) | TWR(8),                   /* Tile  7 */
          TWR(0) | TWR(4),                   /* Tile  6 */
          TWR(0),                            /* Tile  5 */
          
          TWR(C),                            /* Tile  4 */
          TWR(8) | TWR(C),                   /* Tile  3 */
          TWR(4) | TWR(8),                   /* Tile  2 */
          TWR(0) | TWR(4),                   /* Tile  1 */
          TWR(0)                             /* Tile  0 */

       },


       /* --------------------------------------------- */
       /* Shadowed towers by XY tiles                   */
       /* --------------------------------------------- */
       {
          TWR(F),                            /* Tile 24 */
          TWR(E) | TWR(F),                   /* Tile 23 */
          TWR(D) | TWR(E),                   /* Tile 22 */
          TWR(C) | TWR(D),                   /* Tile 21 */
          TWR(C),                            /* Tile 20 */
          
          TWR(B) | TWR(F),                   /* Tile 19 */
          TWR(A) | TWR(B) | TWR(E) | TWR(F), /* Tile 18 */
          TWR(9) | TWR(A) | TWR(D) | TWR(E), /* Tile 17 */
          TWR(8) | TWR(9) | TWR(C) | TWR(D), /* Tile 16 */
          TWR(8) | TWR(C),                   /* Tile 15 */
          
          TWR(7) | TWR(B),                   /* Tile 14 */
          TWR(6) | TWR(7) | TWR(A) | TWR(B), /* Tile 13 */
          TWR(5) | TWR(6) | TWR(9) | TWR(A), /* Tile 12 */
          TWR(4) | TWR(5) | TWR(8) | TWR(9), /* Tile 11 */
          TWR(4) | TWR(8),                   /* Tile 10 */
          
          TWR(3) | TWR(7),                   /* Tile  9 */
          TWR(2) | TWR(3) | TWR(6) | TWR(7), /* Tile  8 */
          TWR(1) | TWR(2) | TWR(5) | TWR(6), /* Tile  7 */
          TWR(0) | TWR(1) | TWR(4) | TWR(5), /* Tile  6 */
          TWR(0) | TWR(4),                   /* Tile  5 */
          
          TWR(3),                            /* Tile  4 */
          TWR(2) | TWR(3),                   /* Tile  3 */
          TWR(1) | TWR(2),                   /* Tile  2 */
          TWR(0) | TWR(1),                   /* Tile  1 */
          TWR(0)                             /* Tile  0 */
         },
   };


   /* 
    | Form the complete set by OR'ing the list from each face 
    | Note that bit mask for the XY plane plane needs to be
    | right justified, that's the '<< 7' stuff.
   */
//   printf("XZ_IGNORE 0x%8.8x  YZ_IGNORE 0x%8.8x\n",~XZ_IGNORE,~YZ_IGNORE);
   return formThrottleByFace (xz & ~XZ_IGNORE, ShadowedTowers.xz)
        | formThrottleByFace (yz & ~YZ_IGNORE, ShadowedTowers.yz)
        | formThrottleByFace (xy << 7,        ShadowedTowers.xy);
}





static unsigned int formThrottleByFace (int                        tiles, 
                                        const unsigned short int *towers)
/*
   DESCRIPTION
   -----------
   Computes the shadowed towers for the specified tiles. For each set
   bit in the 32 bit tile mask, the appropriate element of the 
   shadowedTower array is fetched. This gives the list of towers 
   considered to be shadowed by that tile. The shadowed tower list
   is such that element 0, corresponds to bit 31 being set, that is
   the list is in big endian order.

   PARAMETERS
   ----------
        tiles: A rightt justified bit mask of the struck towers. That
               is bit 31 refers to the highest number tile in the
               face(s) the tile mask is representing. For example,
               for the XY (top) face, tile 24 (numbering from 0)
               occupies bit 31 and tile 0 occupies bit 7.

       towers: An array giving the list (as a bit mask) of the
               shadowed towers by tile. This is read only data.

  RETURNS
  -------
  The list of shadowed tiles, with bit 0 representing tower 0.
*/
{
   int itile    = 0;
   int throttle = 0;


   /*
    |  Loop over the list of struck tiles, 1 per bit.
    | 
    |  !!! IMPROVEMENT !!!
    |  ------------------- 
    |  A better implementation would be to use a FFS instruction (Find First
    |  Set). Almost every modern CPU has such an instruction that locates
    |  the first set bit in an integer in single cyclye. However, I did not
    |  want to venture into portability issues right away.
    */   
//   printf("Tile Word 0x%8.8x\n",tiles);
   while (tiles)
   {
       /* If the bit is hi, OR in the list of shadowed towers for this tile */
       if (tiles < 0) {
            throttle |= towers[itile];
//            printf("For Tile %i the Towers are 0x%8.8x\n",itile,towers[itile]);
       }  
       itile++;
       tiles <<= 1;
   }


   return throttle;
}
