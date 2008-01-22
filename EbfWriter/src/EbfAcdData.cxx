#include <stdio.h>
#include "EbfAcdData.h"
#include "Event/Digi/AcdDigi.h"


/*
 |  History
 |
 |  When    Who  What
 |  ----    ---  ---------------------------
 |  6/05/03 jjr  Removed reference to ObjectVector.h, the AcdDigi class
 |               internally references information in this include, but
 |               did not include it itself.
 |  5/23/03 jjr  Incorrect mapping for tile 430,b and 014,b
*/


static unsigned int *halfswap (unsigned short int *beg, 
                               unsigned short int *end);


/**
 *
 *  @enum   AcdBrds
 *  @brief  Defines the readout order of the 12 ACD electronics boards
 * 
 *   !!! KLUDGE !!! 
 *   --------------
 *   This is purely a guess at this point. The ACD has not specified the
 *   board readout order, or is it GASU designer that hasn't specified 
 *   the readout order. Anyway, it is not known at this time.
 *
 */
enum AcdBrds
{
/* Modified to adhere to the final AEM cable mapping - blw */
    BRD_1LA =  0,  /*!< Acd readout board 1LA readout order =  0 */
    BRD_2LA =  2,  /*!< Acd readout board 1LA readout order =  2 */
    BRD_2RA =  4,  /*!< Acd readout board 2RA readout order =  4 */
    BRD_3LA =  6,  /*!< Acd readout board 3LA readout order =  5 */
    BRD_4LA =  8,  /*!< Acd readout board 4LA readout order =  8 */
    BRD_4RA =  10,  /*!< Acd readout board 4RA readout order = 10 */

    BRD_1RB =  1,  /*!< Acd readout board 1RB readout order =  1 */
    BRD_2LB =  3,  /*!< Acd readout board 1LB readout order =  2 */
    BRD_2RB =  5,  /*!< Acd readout board 2RB readout order =  5 */
    BRD_3RB =  7,  /*!< Acd readout board 3RB readout order =  7 */
    BRD_4LB =  9,  /*!< Acd readout board 4LB readout order =  9 */
    BRD_4RB = 11   /*!< Acd readout board 4RB readout order = 11 */
};




/**
 *
 *  @struct BrdChn
 *  @brief  Mapping for an ADC fiber into an electronics board and 
 *          channel
 *
 *   This is somewhat of a kludge right now since the current version
 *   of GlastSim does not support redundant readout channels. The
 *   kludge is to map a single ACD tile to 2 readout channels. A 
 *   readout channel is defined by an electronics board and channel
 *   number.
 *
 */
struct BrdChn
{
    unsigned int brdA:8; /*!< Electronics board for A fiber */
    unsigned int chnA:8; /*!< Channel number    for A fiber */
    unsigned int brdB:8; /*!< Electronics board for B fiber */
    unsigned int chnB:8; /*!< Channle number    for B fiber */
};




/**
 *
 *  @def   BC(_brdA, _chnA, _brdB, _chnB)
 *  @brief Fills one tile with its mapping to electronics board and 
 *         channel.
 *
 */
#define BC(_brdA, _chnA, _brdB, _chnB) { _brdA, _chnA, _brdB, _chnB }


/**
 *  @def   BC_UNASSIGNED
 *  @brief The unassigned bit are currently set to 0. This macro
 *         just ensures that unassigned bits pick up a 0 by pointing
 *         them at an unused channels. 
 *
 */
#define BC_UNASSIGNED  BC(BRD_1LA, 1, BRD_1RB, 16)


/**
 *
 *  @var   TileToBrdChan
 *  @brief Defines the mapping of the tiles to an electronics board
 *         and channel.
 *
 *   This map used as input the ACD table Rev 2.0, 11/08/02 as reprinted
 *   in the gem Electronics Module. The unconnected input channels were
 *   assigned outputs by MEH in Table 4 of the same document (see p 30).
 *
 */
static const struct BrdChn TileToBrdChn[4][32] =
{ 
  /* XZ+ | XZ- */
  {  BC(BRD_2LA,  5, BRD_2LB, 12),   // XZM__0  200
     BC(BRD_2LA, 11, BRD_2LB,  6),   // XZM__1  201
     BC(BRD_2RA,  3, BRD_2RB, 14),   // XZM__2  202
     BC(BRD_2RA,  7, BRD_2RB, 10),   // XZM__3  203
     BC(BRD_2RA, 12, BRD_2RB,  5),   // XZM__4  204
     BC(BRD_2LA,  3, BRD_2LB, 14),   // XZM__5  210
     BC(BRD_2LA, 10, BRD_2LB,  7),   // XZM__6  211 
     BC(BRD_2RA,  2, BRD_2RB, 15),   // XZM__7  212
     BC(BRD_2RA,  8, BRD_2RB,  9),   // XZM__8  213
     BC(BRD_2RA, 14, BRD_2RB,  3),   // XZM__9  214
     BC(BRD_2LA,  2, BRD_2LB, 15),   // XZM_10  220
     BC(BRD_2LA,  9, BRD_2LB,  8),   // XZM_11  221
     BC(BRD_2RA,  0, BRD_2RB, 17),   // XZM_12  222
     BC(BRD_2RA,  9, BRD_2RB,  8),   // XZM_13  223
     BC(BRD_2RA, 15, BRD_2RB,  2),   // XZM_14  224
     BC(BRD_2RA, 17, BRD_2LB, 17),   // XZM_15  230

     BC(BRD_4RA, 12, BRD_4RB,  5),   // XZP__0  400
     BC(BRD_4RA,  7, BRD_4RB, 10),   // XZP__1  401
     BC(BRD_4RA,  3, BRD_4RB, 14),   // XZP__2  402
     BC(BRD_4LA, 12, BRD_4LB,  5),   // XZP__3  403
     BC(BRD_4LA,  6, BRD_4LB, 11),   // XZP__4  404
     BC(BRD_4RA, 14, BRD_4RB,  3),   // XZP__5  410
     BC(BRD_4RA,  8, BRD_4RB,  9),   // XZP__6  411
     BC(BRD_4RA,  2, BRD_4RB, 15),   // XZP__7  412
     BC(BRD_4LA, 11, BRD_4LB,  6),   // XZP__8  413
     BC(BRD_4LA,  5, BRD_4LB, 12),   // XZP__9  414
     BC(BRD_4RA, 15, BRD_4RB,  2),   // XZP_10  420
     BC(BRD_4RA,  9, BRD_4RB,  8),   // XZP_11  421
     BC(BRD_4RA,  0, BRD_4RB, 17),   // XZP_12  422
     BC(BRD_4LA, 10, BRD_4LB,  7),   // XZP_13  423
     BC(BRD_4LA,  3, BRD_4LB, 14),   // XZP_14  424
     BC(BRD_4RA, 17, BRD_4LB, 17)    // XZP_15  430
  },
  /*  YZ+ | YZ- */
  {
     BC(BRD_2LA,  1, BRD_1RB,  3),   // YZM__0  100
     BC(BRD_1LA,  6, BRD_1RB,  7),   // YZM__1  101
     BC(BRD_1LA,  9, BRD_1RB,  8),   // YZM__2  102
     BC(BRD_1LA, 10, BRD_1RB, 11),   // YZM__3  103
     BC(BRD_1LA, 14, BRD_4RB,  1),   // YZM__4  104
     BC(BRD_2LA,  0, BRD_1RB,  2),   // YZM__5  110
     BC(BRD_1LA,  5, BRD_1RB,  6),   // YZM__6  111
     BC(BRD_1LA,  8, BRD_1RB,  9),   // YZM__7  112
     BC(BRD_1LA, 11, BRD_1RB, 12),   // YZM__8  113
     BC(BRD_1LA, 15, BRD_4RB,  0),   // YZM__9  114
     BC(BRD_1LA,  0, BRD_1RB,  1),   // YZM_10  120
     BC(BRD_1LA,  4, BRD_1RB,  5),   // YZM_11  121
     BC(BRD_1LA,  7, BRD_1RB, 10),   // YZM_12  122
     BC(BRD_1LA, 12, BRD_1RB, 13),   // YZM_13  123
     BC(BRD_1LA, 16, BRD_1RB, 17),   // YZM_14  124
     BC(BRD_1LA, 17, BRD_1RB,  0),   // YZM_15  130

     BC(BRD_3LA, 14, BRD_2RB,  1),   // YZP__0  300
     BC(BRD_3LA, 10, BRD_3RB, 11),   // YZP__1  301
     BC(BRD_3LA,  9, BRD_3RB,  8),   // YZP__2  302
     BC(BRD_3LA,  6, BRD_3RB,  7),   // YZP__3  303
     BC(BRD_4LA,  1, BRD_3RB,  3),   // YZP__4  304
     BC(BRD_3LA, 15, BRD_2RB,  0),   // YZP__5  310
     BC(BRD_3LA, 11, BRD_3RB, 12),   // YZP__6  311
     BC(BRD_3LA,  8, BRD_3RB,  9),   // YZP__7  312
     BC(BRD_3LA,  5, BRD_3RB,  6),   // YZP__8  313
     BC(BRD_4LA,  0, BRD_3RB,  2),   // YZP__9  314
     BC(BRD_3LA, 16, BRD_3RB, 17),   // YZP_10  320
     BC(BRD_3LA, 12, BRD_3RB, 13),   // YZP_11  321
     BC(BRD_3LA,  7, BRD_3RB, 10),   // YZP_12  322
     BC(BRD_3LA,  4, BRD_3RB,  5),   // YZP_13  323
     BC(BRD_3LA,  0, BRD_3RB,  1),   // YZP_14  324
     BC(BRD_3LA, 17, BRD_3RB,  0),   // YZP_15  330
    },
    /* TOP */
  {
     BC(BRD_2LA,  6, BRD_2LB, 11),   // TOP__0 000
     BC(BRD_2LA, 12, BRD_2LB,  5),   // TOP__1 001
     BC(BRD_2LA, 17, BRD_2LB,  0),   // TOP__2 002
     BC(BRD_2RA,  6, BRD_2RB, 11),   // TOP__3 003
     BC(BRD_2RA, 11, BRD_2RB,  6),   // TOP__4 004
     BC(BRD_2LA,  7, BRD_2LB, 10),   // TOP__5 010
     BC(BRD_2LA, 13, BRD_2LB,  4),   // TOP__6 011
     BC(BRD_2RA,  4, BRD_2RB, 13),   // TOP__7 012
     BC(BRD_2RA,  5, BRD_2RB, 12),   // TOP__8 013
     BC(BRD_2RA, 10, BRD_2RB,  7),   // TOP__9 014
     BC(BRD_2LA,  8, BRD_2LB,  9),   // TOP_10 020
     BC(BRD_2LA, 14, BRD_2LB,  3),   // TOP_11 021
     BC(BRD_2LA, 15, BRD_2LB,  2),   // TOP_12 022
     BC(BRD_4LA, 15, BRD_4LB,  2),   // TOP_13 023
     BC(BRD_4LA,  9, BRD_4LB,  8),   // TOP_14 024
     BC(BRD_4RA, 10, BRD_4RB,  7),   // TOP_15 030
     BC(BRD_4RA,  5, BRD_4RB, 12),   // TOP_16 031
     BC(BRD_4RA,  4, BRD_4RB, 13),   // TOP_17 032
     BC(BRD_4LA, 14, BRD_4LB,  3),   // TOP_18 033
     BC(BRD_4LA,  8, BRD_4LB,  9),   // TOP_19 034
     BC(BRD_4RA, 11, BRD_4RB,  6),   // TOP_20 040
     BC(BRD_4RA,  6, BRD_4RB, 11),   // TOP_21 041
     BC(BRD_4LA, 17, BRD_4LB,  0),   // TOP_22 042
     BC(BRD_4LA, 13, BRD_4LB,  4),   // TOP_23 043
     BC(BRD_4LA,  7, BRD_4LB, 10),   // TOP_24 044

     BC_UNASSIGNED               ,   // Unassigned_25
     BC_UNASSIGNED               ,   // Unassigned_26
     BC_UNASSIGNED               ,   // Unassigned_27
     BC_UNASSIGNED               ,   // Unassigned_28
     BC_UNASSIGNED               ,   // Unassigned_29
     BC_UNASSIGNED               ,   // Unassigned_30
     BC_UNASSIGNED               ,   // Unassigned_31
    },

    /* RIBBONS and UNASSIGNED ELECTRONICS CHANNELS */
    {
     BC(BRD_3LA, 13, BRD_1RB,  4),  // RBN_04 500
     BC(BRD_3LA,  2, BRD_1RB, 15),  // RBN_05 501
     BC(BRD_1LA,  2, BRD_3RB, 15),  // RBN_06 502 
     BC(BRD_1LA, 13, BRD_3RB,  4),  // RBN_07 503
     BC(BRD_2LA,  4, BRD_4RB,  4),  // RBN_00 600
     BC(BRD_4RA,  1, BRD_2LB,  1),  // RBN_01 601
     BC(BRD_2RA,  1, BRD_4LB,  1),  // RBN_02 602 
     BC(BRD_4LA,  4, BRD_2RB,  4),  // RBN_03 603 

     BC_UNASSIGNED               ,   // Unassigned_08
     BC_UNASSIGNED               ,   // Unassigned_09
     BC_UNASSIGNED               ,   // Unassigned_10
     BC_UNASSIGNED               ,   // Unassigned_11
     BC_UNASSIGNED               ,   // Unassigned_12
     BC_UNASSIGNED               ,   // Unassigned_13
     BC_UNASSIGNED               ,   // Unassigned_14
     BC_UNASSIGNED               ,   // Unassigned_15

     BC(BRD_4LA,  2, BRD_4LB, 15),   // Unassigned_16 ( 0)
     BC(BRD_4LA, 16, BRD_4LB, 16),   // Unassigned_17 ( 1)
     BC(BRD_4RA, 13, BRD_4RB, 16),   // Unassigned_18 ( 2)
     BC(BRD_4RA, 16, BRD_2LB, 16),   // Unassigned_19 ( 3)
     BC(BRD_1LA,  1, BRD_1RB, 16),   // Unassigned_20 ( 4)
     BC(BRD_1LA,  3, BRD_1RB, 14),   // Unassigned_21 ( 5)
     BC(BRD_2LA, 16, BRD_2LB, 13),   // Unassigned_22 ( 6)
     BC(BRD_2RA, 13, BRD_2RB, 16),   // Unassigned_23 ( 7)

     BC(BRD_3LA,  3, BRD_3RB, 14),   // Unassigned_24 ( 8)
     BC(BRD_3LA,  1, BRD_3RB, 16),   // Unassigned_25 ( 9)
     BC(BRD_2RA, 16, BRD_4LB, 13),   // Unassigned_26 (10)
     BC_UNASSIGNED               ,   // Unassigned_27
     BC_UNASSIGNED               ,   // Unassigned_28
     BC_UNASSIGNED               ,   // Unassigned_29
     BC_UNASSIGNED               ,   // Unassigned_30
     BC_UNASSIGNED               ,   // Unassigned_31
    }
};


    


/**
 *
 *  @fn    initialize ()
 *  @brief Initializes a EbfAcdData to handle a new event
 *
 */
void EbfAcdData::initialize ()
{
    m_gem.cno         = 0;

    m_gem.vetoes [XZ] = 0;
    m_gem.vetoes [YZ] = 0;
    m_gem.vetoes [XY] = 0;
    m_gem.vetoes [RU] = 0;

    m_gem.accepts[XZ] = 0;
    m_gem.accepts[YZ] = 0;
    m_gem.accepts[XY] = 0;
    m_gem.accepts[RU] = 0;

    for (int ibrd = 0; ibrd < NumBoards; ibrd++)
    {
        m_brds[ibrd].accepts = 0;
        m_brds[ibrd].vetoes  = 0;
        for(int ichn = 0; ichn < NumChannelsPerBoard; ichn++)
        {
           m_brds[ibrd].adcs[ichn] = 0;
           m_brds[ibrd].adcsRng[ichn] = 0;
        }
    }
    
}





/**
 *
 *  @fn     fill (const AcdDigiCol &acd)
 *  @brief  Fills the ACD hits and the ADC values.
 *
 */
void EbfAcdData::fill (const Event::AcdDigiCol &tiles)
{
   /* 
    | The hardware packs the ACD channel bits in big endian fashion, that
    | is, channel 0 goes into the most significant bit of the 18 bit mask.
   */
#  define BIT_17 (1<<17)

   /* Initialize the data structure for this event */
   initialize ();


   /* Protect against an empty ACD record, improbable but does occur */
   if (tiles.size ())
   {
       unsigned int       cno_t = 0;
       unsigned int     *vetoes = m_gem.vetoes;
       unsigned int    *accepts = m_gem.accepts;
            
       
       for (Event::AcdDigiCol::const_iterator it = tiles.begin();
            it != tiles.end();
            it++)
       {
           const Event::AcdDigi& digi = **it;
           idents::AcdId           id = digi.getId();
           unsigned int         pha_a = digi.getPulseHeight(Event::AcdDigi::A);
           unsigned int         pha_b = digi.getPulseHeight(Event::AcdDigi::B);
           unsigned int         rng_a = digi.getRange(Event::AcdDigi::A);
           unsigned int         rng_b = digi.getRange(Event::AcdDigi::B);
           bool                veto_a = digi.getVeto       (Event::AcdDigi::A);
           bool                veto_b = digi.getVeto       (Event::AcdDigi::B);
           bool                 cno_a = digi.getHighDiscrim(Event::AcdDigi::A);
           bool                 cno_b = digi.getHighDiscrim(Event::AcdDigi::B);
           bool              accept_a = digi.getLowDiscrim (Event::AcdDigi::A);
           bool              accept_b = digi.getLowDiscrim (Event::AcdDigi::B);
           
/*
           printf("Channel     A             B\n");
           printf("pha      0x%8.8x     0x%8.8x\n",pha_a,pha_b);
           printf("veto     0x%8.8x     0x%8.8x\n",veto_a,veto_b);
           printf("cno      0x%8.8x     0x%8.8x\n",cno_a,cno_b);
           printf("accept   0x%8.8x     0x%8.8x\n",accept_a,accept_b);
           printf("------------------------------\n");
*/
//  cno on side B never seems to be set; so do it to match side A
///           cno_b = cno_a;           


           int             oface;
           struct BrdChn  brdChn;
           struct AcdBrd   *brdA;
           struct AcdBrd   *brdB;
           int              chnA;
           int              chnB;

//           id.Print("C");
           
           int            face = id.face   ();
           int           channel = 0;
           if(id.tile()) {
               channel = (5 * id.row() + id.column());
//               printf("EbfAcd Hit Face %i Row %i  Column %i  Chan %i\n",face,id.row(),id.column(),channel);
           } else if(id.ribbon()) {            
//               printf("EbfAcd Hit Face %i Ribbon Orientation %i Ribbon Number %i\n",face,id.ribbonOrientation(),id.ribbonNum(),channel);
               channel = id.ribbonNum();
           } else {
               // Suppress these warnings after 5 times.  Real data, actually does contain "data" from N/A's
               // we are not particularly interested in this message in that case.
               static int warnCount=0;
               if (warnCount < 5)
                   printf("WARNING: ACD PMT that is not a Tile or Ribbon\n");
               else if (warnCount == 5)
                   printf("WARNING: ACD PMT is not Tile or Ribbon messages will be suppressed from this point\n");
               ++warnCount;
           }

           /*
            | Compute the output face number and bit offset, adjust mask
            | I could not find an enumeration defining the face numbering.
            | This has been gleaned from the ACD tile map document.
            |
            | Not that the 5 and 6 denote the ribbons, with 5 representing
            | the A ribbons and 6 the B ribbons. The 4 A ribbons occupy the
            | low 4 bits and the B ribbons bits the next 4 bits. GLEAM 
            | currently does not generate and values on the unassigned
            | channels. While there are no tiles connected to these channels
            | they readouts are connected to front-end electronics, which,
            | in principle (due to cross-talk for example) could produce
            | a signal.
           */
           if      (face == 0) { oface = XY;                }
           else if (face == 1) { oface = YZ;                }
           else if (face == 2) { oface = XZ;                }
           else if (face == 3) { oface = YZ; channel += 16; }
           else if (face == 4) { oface = XZ; channel += 16; }
           else if (face == 5) { oface = RU;                }
           else if (face == 6) { oface = RU; channel +=  4; }


           brdChn = TileToBrdChn[oface][channel];
           brdA   = m_brds + brdChn.brdA;
           brdB   = m_brds + brdChn.brdB;
           chnA   = brdChn.chnA;
           chnB   = brdChn.chnB;

// Look for entries in the Not Assigned Location
/*            if(oface==RU) {
                printf("Face: %i   oFace: %i \n",face,oface);
                printf("Channel %i    A             B\n",channel);
                printf("pha        0x%8.8x     0x%8.8x\n",pha_a,pha_b);
                printf("veto       0x%8.8x     0x%8.8x\n",veto_a,veto_b);
                printf("cno        0x%8.8x     0x%8.8x\n",cno_a,cno_b);
                printf("accept     0x%8.8x     0x%8.8x\n",accept_a,accept_b);
                printf("------------------------------\n");                                
//            } 
*/


           /*
            | !!! KLUDGE !!!
            | --------------
            | There is a bug in GLEAM which causes the B side VETO and
            | CNO bits not to be set. For testing purposes, I am just
            | going to hardwire some numbers. This should be removed when
            | GLEAM is fixed. These numbers are all WAGS, but basis of 
            | the WAG is
            |
            |      1 MIP  = ~150 counts
            |      VETO   = .4  MIPs =>  60
            |      ACCEPT = .15 MIPs =>  23 
            |      CNO    =  3  MIPs => 450
            |
            | This is good enough for my immediate purposes, but does need
            | to get fixed.
           
           accept_a = pha_a > 23;
           accept_b = pha_b > 23;
           
           veto_a   = pha_a > 60;
           veto_b   = pha_b > 60;

           cno_a    = pha_a > 450;
           cno_b    = pha_b > 450;
           */

           /* Set the veto bits... */
           if (veto_a)
           {
               vetoes[oface] |= (1      << channel);
               brdA->vetoes  |= (BIT_17 >>    chnA);
           }
           

           if (veto_b)
           {
               vetoes[oface] |= (1      << channel);
               brdB->vetoes  |= (BIT_17 >> chnB);
           }
           

           /* Store the pha and range values */
           brdA->adcsRng[chnA] = rng_a;
           brdA->adcs[chnA] = pha_a;

           brdB->adcsRng[chnB] = rng_b;
           brdB->adcs[chnB] = pha_b;           


           /* If above the ADC accept threshold ... */
           if (accept_a)
           {
               /* Add to the bit mask */
               accepts[oface] |= (1      << channel);
               brdA->accepts  |= (BIT_17 >> chnA);

           }
           

           if (accept_b)
           {
               /* Add to the bit mask */
               accepts[oface] |= (1      << channel);
               brdB->accepts  |= (BIT_17 >> chnB);
           }
           


           /* 
            | Compute the CNO trigger
            | Note that since this is a trigger, the packing style is 
            | little-endian, ie Board 0 goes into the least significant bit.
           */
           if (cno_a)
           {
               cno_t |= (1 << brdChn.brdA);
           }

           if (cno_b)
           {
               cno_t |= (1 << brdChn.brdB);      
           }

  
// Print Information about this board
/*           printf("BrdA %i\n",brdChn.brdA);
           printf("Accepts 0x%8.8x\n",brdA->accepts);
           printf("Vetoes  0x%8.8x\n",brdA->vetoes);  
           for(int chnl=0; chnl<18; chnl++){
              printf("Chnl %2i  ADC 0x%8.8x  Rng 0x%8.8x\n",chnl,brdA->adcs[chnl],brdA->adcsRng[chnl]);              
           }
           printf("BrdB %i\n",brdChn.brdB);
           printf("Accepts 0x%8.8x\n",brdB->accepts);
           printf("Vetoes  0x%8.8x\n",brdB->vetoes);  
           for(int chnl=0; chnl<18; chnl++){
              printf("Chnl %2i  ADC 0x%8.8x  Rng 0x%8.8x\n",chnl,brdB->adcs[chnl],brdB->adcsRng[chnl]);              
           }           
*/           
       }

       m_gem.cno = cno_t;
   }
   
   return;
}




/**
 *
 *  @fn         int format (unsigned int *dst)
 *  @brief      Converts the CAL data into Event Builder Format.
 *  @param  dst The destination array.
 *  @return     Pointer to the next available word in the output array.
 *
 *   It is the responsibility of the caller to ensure that enough room
 *  for a maximum record exists in the output array.  
 *
 */
unsigned int *EbfAcdData::format (unsigned int *dst) const
{
   /* 
    | Gives the outline of an ACD output block for 1 FREE board. The PHA
    | values are packed 16 bits for every accept bit, or until the maximum
    | number of PHAs (~4) have been reached.
    |
    | This information is taken from:
    |
    |      The ACD Electronics Module (AEM)
    |      A Primer
    |      LAT-TD-00039-D1
    |      June 7,2002
    |
    |
    |   3 3 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1
    |   1 0 9 8 7 6 5 4 3 2 1 0 9 8 7 6 5 4 3 2 1 0 9 8 7 6 5 4 3 2 1 0
    |  +-------------+-----------------+-------------+-----------------+
    |  |                                                               |
    |  |S V V V V V V V V V V V V V V V V V V A A A A A A A A A A A A A|
    |  |t e e e e e e e e e e e e e e e e e e c c c c c c c c c c c c c|
    |  |a t t t t t t t t t t t t t t t t t t c c c c c c c c c c c c c|
    |  |r o o o o o o o o o o o o o o o o o o c c c c c c c c c c c c c|
    |  |t                                     e e e e e e e e e e e e e|
    |  |                                      p p p p p p p p p p p p p| WORD 0
    |  |                                      t t t t t t t t t t t t t|
    |  |B                                                              |
    |  |i                     1 1 1 1 1 1 1 1                     1 1 1|
    |  |t 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7 8 9 0 1 2|
    |  |                                                               |
    |  +-------------+-----------------+--------------+----------------+
    |  |                               |         PHA READOUT 0         |
    |  |A A A A A P P E C C C C 0 0 0 0 M R V V V V V V V V V V V V V P|
    |  |c c c c c H a n a a a a         o a a a a a a a a a a a a a a a|
    |  |c c c c c A r d b b b b         r n l l l l l l l l l l l l l r|
    |  |e e e e e   i   l l l l         e g u u u u u u u u u u u u u i|
    |  |p p p p p   t o e e e e           e e e e e e e e e e e e e e t| WORD 1
    |  |t t t t t   y f                                               y|
    |  |                N N N N                                        |
    |  |1 1 1 1 1   E C u u u u                                       E|
    |  |3 4 5 6 7   r a m m m m                                       r|
    |  |            r b b b b b                                       r|
    |  |              l e e e e                                        |
    |  |              e r r r r                                        |
    |  |                0 1 2 3
    |  +-------------+-----------------+--------------+----------------+
    |  |         PHA READOUT 1         |       PHA READOUT 2           | WORD 2
    |  |         PHA READOUT 3         |       PHA READOUT 4           | WORD 3
    |  |                     .         |                   .           |
    |  |                     .         |                   .           |
    |  |         PHA READOUT n-1       |       PHA READOUT n           | WORD n
    |  +-------------+-----------------+--------------+----------------+
    |
    | Note that this must be laid out as a series of 16-bit words. This
    | is truly unforturante. This is the only data from the LAT which is
    | not 32-bit oriented. If it were not for the ACD data, a LAT event
    | in EBF format could be processed on a little endian machine by 
    | simple word swapping the whole event (or whole collection of events)
    | without regards to its internal meaning.
   */
#  define START_BIT (1 << 31)
   unsigned short int  *out = (unsigned short int *)dst;  // Output word 
   unsigned short int  *start = out;  
   const struct AcdBrd *brd = m_brds;                     // ADC board info


   /* Loop over each electronics board */
   for (unsigned int ibrd = 0; ibrd < 12; brd++, ibrd++)
   {
       int               accepts = brd->accepts;
       int                vetoes = brd->vetoes;
       const unsigned short *adc = brd->adcs;
       const unsigned short *rng = brd->adcsRng;
       unsigned int start_v0_17_a0_12;
       unsigned int a13_a17_dat_parE_cabNum;
       
       /* Print out boards and Accepts */
//       printf("Board %i  Vetoes 0x%8.8x  Accepts 0x%8.8x\n",ibrd,vetoes,accepts);


       /* Build the first 32 bits */
       start_v0_17_a0_12 = START_BIT | (vetoes  << 13) | (accepts >> 5);
       
       start = out;
       
       /* Stash the first 2 16 bit words of the fixed portion */
       *out++ = (start_v0_17_a0_12 >> 16);
       *out++ = start_v0_17_a0_12 & 0xffff;
       //printf("EbfACD: 1st 16b header %x\n",(start_v0_17_a0_12 >> 16));
       //printf("EbfACD: 2nd 16b header %x\n",(start_v0_17_a0_12 & 0xffff));

       /* build third 16 bit word of the fixed portion*/
       a13_a17_dat_parE_cabNum = (accepts << 11) | (accepts ? (1<<10) : 0) |
                                 (ibrd << 4);
       /* put last cable flag in...KLUDGE: this assumes they are all there*/                          
       if(ibrd==11) a13_a17_dat_parE_cabNum |= 1<<8;
       
       /* Stash the third 16 bit word */
       *out++ = a13_a17_dat_parE_cabNum & 0xffff;
       //printf("EbfACD: 3rd 16b header %x\n",(a13_a17_dat_parE_cabNum & 0xffff));

//       for(int i=0; i<3;i++) {
//         printf("Raw Word %i  0x%4.4x\n",i,*start); 
//         start++;
//       }

       /* Left justify the accepts */
       accepts <<= 14;


       /* Loop over accept mask for this board */
       while (accepts)
       {
           if ((signed int)accepts < 0)
           {
               unsigned short pha;
               unsigned short rngV;
               /* 
                | The 16 bit PHA word is laid out as follows
                |
                |  f e d c b a 9 8 7 6 5 4 3 2 1 0
                |    M R A A A A A A A A A A A A P
                |
                |  Where
                |      P: Parity Error Bit
                |      A: 12 bit ADC value
                |      R: ADC range bit
                |      M: More ADC to follow, if 0, this is the last one
                |
                | !!! KLUDGE !!!
                | --------------
                | This is pretty inelegant coding. One should lay out a
                | bit structure to define this word. Unfortunately the
                | structure definition depends on a combination of the
                | compiler and machine endianness. Depending on the
                | combination, the bit structure may have to be laid out
                | LSB to MSB or MSB to LSB. GLEAM does not have (or I was
                | unable to locate anyone who knew of) a symbol indicating
                | the lay out. So I've kind of hardcoded this in; not good.
               */
               pha       = *adc << 1;
               rngV      = *rng << 13;
               accepts <<= 1;
               
                          
               /* 
                | If out of ADC/pha values on this board. 
                |
                | NOTE
                | ----
                | The original specification called for the FREE boards to
                | limit the number of ADC values it would emit. This was
                | attempt to contain the deadtime. The FREE boards take 
                | about 10-12 usecs to digitize the data + .8 usecs/ADC to
                | transmit the bits. If every ADC was transmitted (18) 
                | this would take > 30usecs, those raising the deadtime
                | bar from 20usecs. However, this was deemed 'too hard' to
                | implement, so one see little breadcrumbs in the data
                | format needed to support this now non-existent feature
                | For example the 'more' bit in the PHA value. The number
                | PHA must not exactly match the number bits in the accept
                | vector, so this bit is unnecessary, but the hardware still
                | emits it.
               */
               if (accepts == 0) 
               {
                   /* Just store the value and call it quits */
                   *out++ = rngV | pha;
                   //printf("EbfACD: pha value(f) %x\n",pha);
                   break;
               }
               else
               {
                   /* Store this value and indicate more to come */
                   *out++ = pha | 0x4000;
                   //printf("EbfACD: pha value %x\n",pha | 0x4000);
               }
           }
           else
           {
               accepts <<= 1;
           }
           
           adc++;
       }
   }


   /* 16-bit swap and return a pointer to the next available 32-bit word */
   return halfswap ((unsigned short int *)dst, out);
}


/**
 *
 *  @fn     void print () const
 *  @brief  Diagnostic print routine 
 * 
 *  Routine to print the hits on each ACD tile information. This
 *  is mainly used to debug the filling routine
 *
 */
void EbfAcdData::print() const
{
   std::printf ("\n"
                "ACD data\n"
                "  XZ Tiles = %8.8x %8.8x\n"
                "  YZ Tiles = %8.8x %8.8x\n"
                "  XY Tiles = %8.8x %8.8x\n"
                "  RU Tiles = %8.8x %8.8x\n"
                "\n",
                m_gem.vetoes[XZ], m_gem.accepts[XZ],
                m_gem.vetoes[YZ], m_gem.accepts[YZ],
                m_gem.vetoes[XY], m_gem.accepts[XY],
                m_gem.vetoes[RU], m_gem.accepts[RU]);


   for (unsigned int list = 0; 
        list < sizeof(m_gem.vetoes)/sizeof(*m_gem.vetoes);
        list++)
   {
       static const char FaceNames[8][4] =
           {  "XY ", "YZ-", "XZ-", "YZ+", "XZ+", "RB ", "RA ", "UA "  };
       
       unsigned int      veto = m_gem.vetoes[list];
       unsigned int    accept = m_gem.accepts[list];
       unsigned int      tile = veto | accept;
       unsigned int   channel =             0;

       
       /* Loop over the struck tiles (either veto or accept) on this face */
       do
       {
           int             face;
           char      vetoSymbol;
           char    acceptSymbol;
           struct BrdChn brdChn;
           int             brdA;
           int             brdB;
           int             chnA;
           int             chnB;
           int              row;
           int              col;
           

           /* Is either this accept or veto tile struck */
           if ( (tile & 1) == 0) continue;

           /* Use a symbol to indicate which (accept or veto) is struck */
           vetoSymbol   = (veto   & 1) ? 'V' : ' ';
           acceptSymbol = (accept & 1) ? 'A' : ' ';
           brdChn       = TileToBrdChn[list][channel];
           brdA         = brdChn.brdA;
           brdB         = brdChn.brdB;
           chnA         = brdChn.chnA;
           chnB         = brdChn.chnB;


           if (list != RU)
           {
               /* Map to standard ACD LAT face numbers */
               static const char Face[6] = { 2, 4, 1, 3, 0, 0};
               int chn;
               
               chn   = (list != XY && channel>=16) ? channel-16 : channel;
               face  = Face[2*list + channel/16];
               row   = chn / 5;
               col   = chn % 5;
           }
           else
           {
               /*
                | Upper 16 bits of RU = Unassigned
                |      bits 4-7 of RU = RA (ribbon A's)
                |      bits 0-3 of RU = RB (ribbon B's)
               */
               col = channel >= 16 ?  face = 7, channel - 16 // Unassigned
                   : channel >=  4 ?  face = 6, channel - 4  // RA's
                                   :  face=  5, channel;     // RB's
               row = 0;
           }
           
           std::printf (" %s[%1d:%1d:%1d] :"
                        "  Brd[%2d:%2d] = %4.4x"
                        "  Brd[%2d:%2d] = %4.4x %c %c\n",
                        FaceNames[face], face, row, col,
                        brdA, chnA, m_brds[brdA].adcs[chnA],
                        brdB, chnB, m_brds[brdB].adcs[chnB],
                        acceptSymbol, vetoSymbol);
       }
       while (channel++, veto >>= 1, accept >>=1, tile >>= 1);

   }
   
   return;
}




static unsigned int *halfswap (unsigned short int *beg, 
                               unsigned short int *end)
/*
   DESCRIPTION
   -----------
   Swaps all the 16 bit words in the specified range. 

   PARAMETERS
   ----------
          beg: Pointer to the first 16 bit word to swap. 
               This must be 32 bit aligned.

          end: Pointer to the last  16 bit word to swap
               (actually 1 past it.) This may or may not
               be 32 bit aligned.

   RETURNS
   -------
   Pointer to the next available 32 bit word 
*/
{
    int                 len =  end - beg;
    unsigned short int *cur = beg;
    

    /* If the length is odd, pad to 32 bit boundary */
    if (len & 1)
    {
        /* Fill in a fun word for the padding */
        *end = 0x0;
        len++;
    }


    len /= 2;
    while (--len >= 0)
    {
        unsigned short int tmp0 = cur[0];
        unsigned short int tmp1 = cur[1];
        
        cur[0] = tmp1;
        cur[1] = tmp0;
        
        cur += 2;
    }
    
    return (unsigned int *)cur;
}



