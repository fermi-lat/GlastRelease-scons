// $Header$
#ifndef CalibData_BadStrips_h
#define CalibData_BadStrips_h

#include "CalibData/CalibBase.h"
#include <vector>

/**
         @file CalibData/BadStrips.h

         This file actually defines two public classes in the CalibData
         namespace : BadStrips and BadStripsVisitor.  BadStrips is
         the TCDS form of bad strip (dead or hot but not both) calibration
         information.  BadStripsVisitor is a pure virtual class to be 
         used as a base by applications wishing to access the information 
         in BadStrips.

         @author  J. Bogart
*/
namespace CalibData {

  /** Visitor callbacks can indicate whether traversal should continue
       or not.
         CONT        normal return: continue processing
         USER_DONE   client has all information desired; no more traversal
         ERROR       client has serious error; abort
         DONE        not used by client.  Will be returned by 
                     BadStrips::traverse   in case processing was normal.
  */
  enum eVisitorRet {CONT, USER_DONE, ERROR, DONE};

  typedef std::vector<unsigned short int> StripCol;

  /**
         Concrete implementation of BadStripsVisitor is required to
         traverse TCDS BadStrips.  Two functions are required:  one
         for (entirely bad) towers and another for (entirely or partially)
         bad uniplanes.  The visitor may be called back more than once
         for the same plane if there are strips with different "degrees of
         badness" in the plane.
  */
  class BadStripsVisitor {
  public:
      
    /**  Handle bad tower
        @param row         zero-based row of tower
        @param col         zero-based column of tower
        @param badness     if gradations of badness have been recorded,
                           larger number here corresponds to worse failure
    */
    virtual eVisitorRet badTower(unsigned int row, unsigned int col, 
                                 int badness)=0;

    /**  Handle bad uniplane with some or all bad strips
        @param row         zero-based row of tower
        @param col         zero-based column of tower
        @param badness     if gradations of badness have been recorded,
                           larger number here corresponds to worse failure
        @param strips      vector of strips of badness @arg badness.  If
                           empty, entire plane is bad.
    */
    virtual eVisitorRet badPlane(unsigned int row, unsigned int col, 
                                 unsigned int tray, bool top,
                                 int badness, 
                                 const StripCol& strips)=0;

  };    // end pure virtual visitor class definition

/** @class BadStrips
    @brief TCDS representation of a collection of hot or dead strips
           (but not both).
*/
  class BadStrips :  public CalibBase {
    /// Allow converter access to private build, set stuff
    friend class XmlBadStripsCnv;  
    
  public:
    enum eBadType {
      DEAD = 0,
      HOT  = 1 };

    eBadType getBadType() const;

    /// call back method for client to access all data
    eVisitorRet traverse(BadStripsVisitor *visitor) const;

    /// Make a contentless (except for CalibBase stuff and hot/dead
    /// discriminator) class.
    BadStrips(eBadType bType, const ITime& since, const ITime& til, 
              int serNo = -1);
    virtual ~BadStrips();

    // Re-implemented from DataObject
    inline virtual const CLID& clID() const { return classID(); }
    
    static inline const CLID& classID() { return CLID_Calib_TKR_BadChan; }
  private:
    class Uniplane {
    public:
      Uniplane(bool allBad, int howBad, int tray, bool top,
               const StripCol& strips);

      Uniplane(const Uniplane& other);

      ~Uniplane();

      /// Short-hand way to declare all strips bad
      bool m_allBad;
      /// Following field is > 0
      int  m_howBad;    
      /// Tray number; zero-based
      unsigned int  m_tray;
      /// m_top is true for top, false for bottom
      bool m_top;
      StripCol*    m_badStrips;
    };

    class Tower {
    public:
      Tower(bool allBad, int howBad, unsigned row, unsigned col);

      Tower(const Tower& other);

      ~Tower();

      /// Short-hand way to declare all strips within the tower bad
      bool m_allBad;
      /// If m_allBad is true, following field is > 0
      int  m_howBad;    
      /// Row number of tower in the array; zero-based
      unsigned short int m_row;
      /// column number of tower in the array; zero-based
      unsigned short m_col;

      std::vector<Uniplane>* m_uniplanes;
    };


    // Re-implemented from CalibBase
    virtual void    update(CalibBase& other);

    // Following are intended for converters to call
    StatusCode addBadTower(bool allbad, int howBad, 
                           unsigned row, unsigned col);

    StatusCode addBadPlane(unsigned short row, unsigned short col,
                           unsigned int tray, bool top, int howBad,
                           StripCol& badStrips);

    // internal stuff
    Tower* findTower(unsigned row, unsigned col);

    eBadType m_type;
    std::vector<Tower>* m_towers;
  };     // end Badstrips class definition
}        // end CalibData namespace
#endif
