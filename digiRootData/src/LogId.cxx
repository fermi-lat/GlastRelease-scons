///////////////////////////////////////////////////////////////////////////
//                                                                        
// The LogId class contains a log ID number (encoding of detector
// element position) as well as functions to obtain the location of
// log in question. LogId is used in the CalDigi class.
//                                                                        
///////////////////////////////////////////////////////////////////////////

// access and set functions for binary calorimeter tags
// Version 1.0 21 Oct 1998 Richard creation
// Version 1.1 25 Oct 1999 R.Dubois Clone from LCD

#include "digiRootData/LogId.h"

ClassImp(LogId)

LogId::LogId(){
    // Default constructor
    // All fields in tagword initialized to 0
    setTag(static_cast<UInt_t>(0));
}

LogId::LogId(UInt_t tag=0) {
    // Create a LogId object with number tag
    setTag(tag);
}

LogId::LogId(UShort_t tower, UShort_t layer, UShort_t column) {
    // Create a log id given the layer, column and tower
    setTag(static_cast<UInt_t>(0));
    setLayer(layer);
    setColumn(column);
    setTower(tower);
}

LogId::~LogId() {
    // Destructor
}


Bool_t LogId::isValidColumn(UInt_t columnVal) {
    // Check if given value is a valid Column number:
    return (
            ((columnVal <= BOUNDS_COLUMN) && !(columnVal & ~CAL_M_COLUMN)) ?
            kTRUE : kFALSE
           );
}

Bool_t LogId::isValidLayer(UInt_t layerVal) {
    // Check if given value is a valid Layer number:
    return (
            ((layerVal <= BOUNDS_LAYER) && !(layerVal & ~CAL_M_LAYER)) ?
            kTRUE : kFALSE
           );
}


Bool_t LogId::isValidTagStruct(TAG_STRUCT ts) {
    // Test the validity of a particular tag structure.
    // Ensures that fields which don't use all possible values are
    // within bounds. 
    if (!isValidColumn(ts.column)) return kFALSE;
    if (!isValidLayer(ts.layer)) return kFALSE;
    
    // If we got this far, everything must be ok...
    return kTRUE;
}

Bool_t LogId::isValidTagWord(UInt_t tagWord) {
    
    TAG_STRUCT ts;
    fillTagStruct(tagWord, &ts);
    
    return isValidTagStruct(ts);
}

// Fill a TAG_STRUCT with data from tagVal.
// *** NOTE: There is no assertion that user supplied tagword contains
// valid fields!
void LogId::fillTagStruct(UInt_t tagVal, LogId::TAG_STRUCT *ts) {
    ts->column = getColumn(tagVal);
    ts->layer = getLayer(tagVal);
    ts->pipeline = 0;//getPipeline(tagVal);
    ts->sequence = 0;//getSequence(tagVal);
    ts->tower = getTower(tagVal);
    ts->xy = getXY(tagVal);
}

// Fill a tagword from a tag structure:
// *** NOTE: There is no assertion that user supplied tag struct contains
// valid fields!
UInt_t LogId::fillTagWord(LogId::TAG_STRUCT *ts) {
    UInt_t tagWord = 0;
    tagWord += (ts->xy & CAL_M_XY) << CAL_V_XY;
    //tagWord += (ts->sequence & CAL_M_SEQ) << CAL_V_SEQ;
    //tagWord += (ts->pipeline & CAL_M_PIPE) << CAL_V_PIPE;
    // When setting layer bits, shed common XY bit from layer value...
    tagWord += ((ts->layer >> CAL_K_XY) & CAL_M_LAYER) << CAL_V_LAYER;
    tagWord += (ts->column & CAL_M_COLUMN) << CAL_V_COLUMN;
    tagWord += (ts->tower & CAL_M_TOWER) << CAL_V_TOWER;

    return tagWord;
}

Bool_t LogId::setTag(UInt_t tagVal) { 
    // Set value of tag word to tagVal if taVal is valid

    if (isValidTagWord(tagVal)) {
        m_tag = tagVal;
        return kTRUE;
    }
    else return kFALSE;
}

Bool_t LogId::setTag(TAG_STRUCT *ts) {
    // Fill tag word from values in a tag struct.

    if (isValidTagStruct(*ts)) {
        m_tag = fillTagWord(ts);
        return kTRUE;
    }
    else return kFALSE;
}

Bool_t LogId::setTower(UInt_t towerVal) {
    // clear the tower bits
    m_tag &= ~(CAL_M_TOWER << CAL_V_TOWER);
    // set the tower bits to the new value
    m_tag |= towerVal << CAL_V_TOWER;
    return kTRUE;
}

Bool_t LogId::setColumn(UInt_t columnVal) {
    // Sets the column field in the tag word

    if (!isValidColumn(columnVal)) 
        return kFALSE;

    // clear the column bits
    m_tag &= ~(CAL_M_COLUMN << CAL_V_COLUMN);
    // set the column bits to the new value
    m_tag |= columnVal << CAL_V_COLUMN;

    return kTRUE;
}

Bool_t LogId::setXY(CALAxes xyVal) {
    // Sets the XY field in the tag word
    // clear the XY bit
    m_tag &= ~(CAL_M_XY << CAL_V_XY);
    // set the XY bit to the new value
    m_tag |= xyVal << CAL_V_XY;

    return kTRUE;
}

Bool_t LogId::setLayer(UInt_t layerVal) {
    // Sets the layer field in the tag word

    // *** See note in header (in the private enum of tagword fields) 
    // *** as to why the following is done:
    // Strip off LSB, as it goes in the XY field...
    UInt_t xyVal = (layerVal & CAL_M_XY);
    // Now get the layer value minus the XY bit...
    UInt_t lVal = (layerVal >> CAL_K_XY);

    // Check if the remaining bits represent a valid layer value:
    if (lVal & ~CAL_M_LAYER)
        return kFALSE;

    // Set LAYER and XY fields...
    // Clear the layer bits
    m_tag &= ~(CAL_M_LAYER << CAL_V_LAYER);
    // set the layer bits to the new value
    m_tag |= lVal << CAL_V_LAYER;
    // clear the XY bit
    m_tag &= ~(CAL_M_XY << CAL_V_XY);
    // set the XY bit to the new value
    m_tag |= xyVal << CAL_V_XY;

    return kTRUE;
}

Bool_t LogId::setId(UShort_t tower, UShort_t layer, UShort_t column) {
    
    if (setLayer(layer) == kFALSE) return kFALSE;
    if (setColumn(column) == kFALSE) return kFALSE;
    if (setTower(tower) == kFALSE) return kFALSE;

    return kTRUE;
}


UInt_t LogId::getTag() const { 
    // Returns the log number - the whole tag word
    return m_tag;
}

UInt_t LogId::getId(UInt_t tagWord) {
    // Returns the id of this log.  Each log in the calorimeter has a unique
    // id value based upon it's position in the cal.
    // Now using Layer,Column,Tower to determine unique id
   // return static_cast<UInt_t>((tagWord >> CAL_V_ID) & CAL_M_ID);
    return static_cast<UInt_t>((tagWord >> CAL_V_NEWID) & CAL_M_NEWID);
}

UInt_t LogId::getId() const {
    // Returns the unique id of this log based on Tower, Layer, Column
    return getId(m_tag);
}


UShort_t LogId::getTower(UInt_t tagWord) {
    // Returns the GLAST tower this log belongs to.  In the '99/2000
    // beamtest, and for the 2001 Balloon flight this value is always 0.
    return static_cast<UShort_t>((tagWord >> CAL_V_TOWER) & CAL_M_TOWER);
}

UShort_t LogId::getTower() const {
    return getTower(m_tag);
}

UShort_t LogId::getColumn(UInt_t tagWord) {
    // Returns the column of the log in the calorimeter
    // Columns are numberd 0-9.  
    // Columns are numbers with increasing coordinate.
    // If this is an X log (oriented along the X axis), the columns are numbered with 
    // increasing Y coordinate.
    // Similarly for the Y logs (oriented along the Y axis), the columns are numbered with
    // increasing X coordinate.
    return static_cast<UShort_t>((tagWord >> CAL_V_COLUMN) & CAL_M_COLUMN);
}

UShort_t LogId::getColumn() const {
    // return column value from tag word:
    return getColumn(m_tag);
}

UShort_t LogId::getLayer(UInt_t tagWord) {
    // Returns the layer this log is in.  Layers are numbered from
    // 0: (closest to TKR) to 7: (furthest from TKR).
    UInt_t layerVal = (tagWord >> CAL_V_LAYER) & CAL_M_LAYER;
    layerVal <<= CAL_K_XY;
    layerVal |= getXY(tagWord);

    return UShort_t(layerVal);
}

UShort_t LogId::getLayer() const {
    // return layer value from tag word:
    return getLayer(m_tag);
}


UShort_t LogId::getXY(UInt_t tagWord) {
    // Returns the X / Y orientation of this log
    return static_cast<UShort_t>((tagWord >> CAL_V_XY) & CAL_M_XY);
}

LogId::CALAxes LogId::getXY() const {
    // Returns whether this log is in an X log or a Y log
    return (CALAxes)getXY(m_tag);
}


