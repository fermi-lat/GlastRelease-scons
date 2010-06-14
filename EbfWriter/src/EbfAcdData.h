#ifndef EBF_ACD_DATA_H
#define EBF_ACD_DATA_H



//
// At one time the arguments were passed by pointer,
// but now they are passed by reference. This means
// that one can no longer forward reference AcdDigiCol
//
#include "Event/Digi/AcdDigi.h"




/**
 *  @class EbfAcdData Tiles
 *  @brief Summarizes the ACD data for use in producing both the EBF
 *         ACD data and the ACD data used and produced by the GLT.
 *
 *   The ACD electronics services 216 channels, 18 channels on each
 *   of 12 so-called FREE boards. The ACD consists of 97 active
 *   elements, 89 tiles and 8 ribbons. Each active element is
 *   readout twice, for a 2-fold redundance.
 *
 *   The ACD data appears in two places. The first is in the ACD
 *   data itself. This is produced by the AEM boards. It has two
 *   components, a bit mask of the ACD elements that are above
 *   threshold and a data structure giving the pulse heights of
 *   those elements that are above a zero-suppression threshold.
 *   The elements in the ACD data are arranged by how they appear
 *   on the FREE boards. Because of physical constraints, the order
 *   of the channels on the FREE boards is rather scrambled, that
 *   is, there is no correlation with the readout order and the
 *   tile geometry.
 *
 *   The second place the ACD data occurs is in the GLT data. Here
 *   the ACD data is reduced to only the veto bit masks. Furthermore
 *   the 2-bold redundance is collapsed by ORing the A and B channels.
 *   The ACD data is, however, reordered to follow the geometry.
 *   The GLT ACD data consists of 4 32 bit words. The words represent
 *   the 32 channels for the tiles measuring in the XZ, the YZ, the
 *   XY (sometimes called the TOP) planes and the ribbons. The ACD in
 *   the GLT also captures the 12 CNO bits (one per FREE board).
 *
 */
class EbfAcdData
{
  public:

    /**
     * @enum  TileFaces
     * @brief Gives symbolic names to the 4 32 bit words of the ACD data
     *        captured in the GLT data structure.
     *
     *  The symbolic names are chosen to name the dimensions that the
     *  tiles measure. For example, the top plane of ACD tiles is named
     *  XY. The XZ and YZ planes can be thought of as 2 16 bit masks,
     *  one set of 16 bits representing the 'minus' side tiles and the
     *  other set representing the 'plus' sides. In other places these
     *  tiles are referred to as YM, YP, XM, YP and TOP. This names
     *  the axis perpendicular to the plane. The problem with this
     *  notation is that it just reads weird in the code. For example,
     *  then intersecting a track with the tiles in the XZ plane, one
     *  must use the tiles named YM (or YP).
     *
     *  Unfortunately, the last word, capturing the hit masks for the
     *  ribbons and unused channels (RU), does not really correspond
     *  to a plane or a 'face'. Oh well.
     */
    enum TileFaces
    {
      XZ = 0,  /*!< Tiles measuring X and Z dimensions                  */
      YZ = 1,  /*!< Tiles measuring Y and Z dimensions                  */
      XY = 2,  /*!< Tiles measuring X and Y dimensions, the TOP tiles   */
      RU = 3,  /*!< The 8 ribbons + the unused channels                 */
      NumMasks /*!< The number of 32 bit words needed to hold the masks */
    };

    static const int NumBoards           = 12; /*!< 12 electronics boards    */
    static const int NumChannelsPerBoard = 18; /*!< 18 channels/board        */
    static const int NumChannels         = NumChannelsPerBoard * NumBoards;
                                               /*!< Total number of channels */
    
    void          initialize ();
    void          fill       (const Event::AcdDigiCol &tiles);
    unsigned int *format     (unsigned int *dst)      const;
    void          print      ()                       const;

    inline unsigned int  cno           () const { return m_glt.cno;         }

    inline const
           unsigned int *vetoes        () const { return m_glt.vetoes;      }
    inline unsigned int  vetoesXY      () const { return m_glt.vetoes[XY];  }
    inline unsigned int  vetoesXZ      () const { return m_glt.vetoes[XZ];  }
    inline unsigned int  vetoesYZ      () const { return m_glt.vetoes[YZ];  }
    inline unsigned int  vetoesRU      () const { return m_glt.vetoes[RU];  }

    inline const
           unsigned int *accepts       () const { return m_glt.accepts;     }
    inline unsigned int  acceptsXY     () const { return m_glt.accepts[XY]; }
    inline unsigned int  acceptsXZ     () const { return m_glt.accepts[XZ]; }
    inline unsigned int  acceptsYZ     () const { return m_glt.accepts[YZ]; }
    inline unsigned int  acceptsRU     () const { return m_glt.accepts[RU]; }
    
  private:

    /**
     *  
     *   @struct GltData
     *   @brief  The data relevant to producing the ACD's contribution
     *           to the GLT data
     *
     */
    struct GltData
    {
        unsigned int                 cno; /*!< Mask of the 12 CNO triggers */
        unsigned int    vetoes[NumMasks]; /*!< XZ, YZ, XY, RU     vetoes   */
        unsigned int   accepts[NumMasks]; /*!< XZ, YZ, XY, RU     accepts  */
    } m_glt;
    

    /**
     *  
     *   @struct AcdBrd
     *   @brief  The data relevant to producing the ACD data contribution
     *
     */
    struct AcdBrd
    {
        unsigned int    vetoes;             /*!< Veto   map for this board */
        unsigned int   accepts;             /*!< Accept map for this board */
        unsigned short
                 adcs[NumChannelsPerBoard]; /*!<       ADCs for this board */
    } m_brds[NumBoards];
};


#endif




