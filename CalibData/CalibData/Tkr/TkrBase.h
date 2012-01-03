// $Header$
/// @file TkrBase.h
/// @author J. Bogart

#ifndef CalibData_TkrBase_h
#define CalibData_TkrBase_h

#include <vector>
#include "CalibData/CalibBase.h"
// #include "CalibData/Tkr/TkrFinder.h"
#include "idents/TkrId.h"

#define TKRBASE_MAXROW 4
#define TKRBASE_MAXCOL 4
#define TKRBASE_MAXTOWER (TKRBASE_MAXROW * TKRBASE_MAXCOL)

class RootTkrBaseCnv;

namespace CalibData {
  class UniBase;
  class TkrFinder;
  /**
     Each derived, completely implemented tkr calibration class should
     register a suitable factory object which will produce the right
     kind of object of a class derived from UniBase
  */
  class UniFactoryBase {
  public: 
    UniFactoryBase() {}
    virtual ~UniFactoryBase() { };
    virtual UniBase* makeUni();
  };

  /**
     Base class for Tkr calibrations other than bad strips, which have
     a fixed amount of data per Si layer, typically involving a
     per-strip or per-gtfe structure
  */
  class TkrBase : public CalibBase {
    friend class ::RootTkrBaseCnv;

  public:
    /**
       Constructor configures its TkrFinder
    */
    TkrBase(unsigned nTowerRow=4, unsigned nTowerCol=4, unsigned nTray=19,
             bool indirect=true);

    virtual ~TkrBase();
    
    virtual UniBase* getUni(const idents::TkrId& id);

    virtual UniBase* getUni(unsigned towerRow, unsigned towerCol,
                                  unsigned tray, bool top);
    
    virtual bool putUni(UniBase* data, const idents::TkrId& id);

    virtual bool putUni(UniBase* data, unsigned towerRow, 
                        unsigned towerCol, unsigned tray, bool top);

    virtual const std::string* getHwserial(unsigned towerRow, 
                                           unsigned TowerCol)
      const;


    // Get dimensioning information; needed when transforming to 
    // permanent storage...might not need this
    /// Get # tower rows
    unsigned getNTowerRow() const;

    /// Get # tower columns
    unsigned getNTowerCol() const;

    /// Get #  trays
    unsigned getNUnilayer() const;

    /// Get pointer to vector of uni for specified tower.
    std::vector<UniBase*> *getUnis(int iTow) {
      // or maybe RootTkrBaseCnv should provide this service
      return &(m_towers[iTow]->m_unis); 
    }

    /// Get # fe chips / unilayer
    // unsigned getNChip() const {return m_finder->getNChip();}
    virtual const CLID& clID() const = 0;     // must be overridden
    static const CLID& classID();   // shouldn't get called

    // Maybe won't need to be virtual after all
    virtual StatusCode update(CalibBase& other, MsgStream* log);

  protected:
    /**
       @class TkrTower
       Represents a tower's worth of tracker calibration data.
       Identify the tower, keep remaining information by uniplane.
       Different derived classes will have their own class
       for per-uniplane data, derived from UniBase.
    */
    class TkrTower {
      friend class ::RootTkrBaseCnv;
    public:
      unsigned m_iRow;
      unsigned m_iCol;
      std::string m_hwserial;
      std::vector <UniBase* > m_unis;
      TkrTower(unsigned row=0, unsigned col=0, 
               std::string hwserial="") 
      : m_iRow(row), m_iCol(col), m_hwserial(hwserial) 
      {        }
      void resize(unsigned n);
      ~TkrTower();
    };                        // end def TkrTower


    // Array below is only of pointers to towers.  Allocate a tower
    // only when needed.
    TkrTower* makeTower(unsigned iTow, unsigned nUni=38);
      


    TkrFinder* m_finder;
    UniFactoryBase* m_factory;
   

    TkrTower* m_towers[TKRBASE_MAXTOWER];
    /// Default is true: keep vector of pointers to data.  
    /// Else derived class keeps vector of data values and must 
    /// do its own fetching and putting.
    bool m_indirect;

  private:
    static const CLID noCLID;

    /** Due to bug in gcc, gdb can't find symbols in constructors.  
        This method is called by the constructor and does most 
        of the work
    */
    void cGuts(unsigned nTowerRow, unsigned nTowerCol, 
               unsigned nTray);
  };


  
}  
#endif

