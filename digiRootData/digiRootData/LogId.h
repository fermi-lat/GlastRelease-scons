

#ifndef LogId_H
#define LogId_H

#include "TObject.h"

/*! \class LogId
\brief
 access and set functions for calorimeter log ids
 Jun 2001 Heather Kelly renamed LogID to LogId
 Version 1.0 21 Oct 1998 Richard creation
 Version 1.1 25 Oct 1999 R.Dubois Clone from LCD towerID
*/

class LogId : public TObject {
public: 
    typedef struct TAG_STRUCT {
        UInt_t xy, 
               sequence, 
               pipeline, 
               layer, 
               column, 
               tower;
    };
private:
    /*
    LogId tag word layout
    (defined in following enum):
     ________________________________________________________________
    |31| ... |17||16|  |  |13|12|  |  |09|08|07||06|  |04|03|  |01|00|
    |__|_____|__||__|__|__|__|__|__|__|__|__|__||__|__|__|__|__|__|__|
    |           ||   TOWER   |  COLUMN   |LAYER||PIPELINE|SEQUENCE|XY|
    |RESERVERD 4||___________|___________|_____||________|________|__| 
    |FUTURE USE ||         LOG POSITION        ||   LOG ID NUMBER    |
    |___________||_____________________________||____________________| 

    ***** NOTES: *****

      (1) LAYER Field:
      The layer field contains the most significant two bits necessary
      to encode the 8 physical cal layers.  When the 'XY' bit is appended
      to this field as the LSB an actual physical layer number is obtained.
      This is done to prevent the dual storage of data which could result 
      in errors if they were to differ.
    */
    enum {
        // X/Y bit:
        CAL_K_XY = 1,
        CAL_V_XY = 0,
        CAL_M_XY = ((1 << CAL_K_XY) - 1),
        // ADC Sequence #:
        CAL_K_SEQ = 3,
        CAL_V_SEQ = (CAL_V_XY + CAL_K_XY),
        CAL_M_SEQ = ((1 << CAL_K_SEQ) - 1),
        // Pipeline #:
        CAL_K_PIPE = 3,
        CAL_V_PIPE = (CAL_V_SEQ + CAL_K_SEQ),
        CAL_M_PIPE = ((1 << CAL_K_PIPE) - 1),
        // Log ID Number as a whole:
        CAL_K_ID = (CAL_K_XY + CAL_K_SEQ + CAL_K_PIPE),
        CAL_V_ID = CAL_V_XY,
        CAL_M_ID = ((1 << CAL_K_ID) - 1),

        // Log Position (Layer):
        CAL_K_LAYER = 2,
        CAL_V_LAYER = (CAL_V_PIPE + CAL_K_PIPE),
        CAL_M_LAYER = ((1 << CAL_K_LAYER) - 1),
        // Log Position (Column):
        CAL_K_COLUMN = 4,
        CAL_V_COLUMN = (CAL_V_LAYER + CAL_K_LAYER),
        CAL_M_COLUMN = ((1 << CAL_K_COLUMN) - 1),
        // Log Position (Tower):
        CAL_K_TOWER = 4,
        CAL_V_TOWER = (CAL_V_COLUMN + CAL_K_COLUMN),
        CAL_M_TOWER = ((1 << CAL_K_TOWER) - 1),

        // Upper Bounds values for bounded fields:
        // NOTE: no bounds limits for fields which use all possible values.
        BOUNDS_COLUMN = 9,
        BOUNDS_PIPE = 4
    };
    UInt_t m_tag;   // Packed word containing log data

    static Bool_t isValidTagStruct(LogId::TAG_STRUCT ts);
    static Bool_t isValidTagWord(UInt_t tagWord);
    static Bool_t isValidColumn(UInt_t columnVal);
    static Bool_t isValidPipeLine(UInt_t pipeVal);

    static UInt_t getTower(UInt_t tagWord);
    static UInt_t getColumn(UInt_t tagWord);
    static UInt_t getLayer(UInt_t tagWord);

    static UInt_t getId(UInt_t tagWord);
    static UInt_t getPipeline(UInt_t tagWord);
    static UInt_t getSequence(UInt_t tagWord);
    static UInt_t getXY(UInt_t tagWord);
public:
    // orientation of the layer
    // X logs are oriented along the X axis
    // Y logs are oriented along the Y axis
    typedef enum {
        X = 0,
        Y
    } CALAxes;

    LogId();            
    LogId(UInt_t tag);
    virtual ~LogId();

    static Bool_t fillIdFromGeom(LogId::TAG_STRUCT *ts);
    static Bool_t fillGeomFromId(LogId::TAG_STRUCT *ts);
    static UInt_t fillTagWord(LogId::TAG_STRUCT *ts);
    static void fillTagStruct(UInt_t tagVal, LogId::TAG_STRUCT *ts);

    UInt_t getTag()         const;
    UInt_t getColumn()      const;
    UInt_t getLayer()       const;

    UInt_t getId()          const;
    CALAxes getXY()         const;
    UInt_t getPipeline()    const;
    UInt_t getSequence()    const;

    Bool_t setTag(UInt_t tagVal);
    Bool_t setTag(TAG_STRUCT *ts);
    Bool_t setColumn(UInt_t columnVal);
    Bool_t setLayer(UInt_t layerVal);
    Bool_t setXY(CALAxes xyVal);
//    Bool_t setID(UInt_t idVal);

ClassDef(LogId,3)       // CsI log identification and readout information
};
#endif
