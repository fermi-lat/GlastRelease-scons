/// @file TkrTruncationInfo.h

#ifndef __TkrTruncationInfo_H__
#define __TkrTruncationInfo_H__

#include "Event/Recon/TkrRecon/TkrTruncatedPlane.h"

#include "GaudiKernel/DataObject.h"

#include <map>

static const CLID& CLID_TkrTruncationInfo = InterfaceID("TkrTruncationInfo",  0, 0);


namespace Event {
    /// small class to hold map key
    class SortId {
        enum {
            shiftVIEW  = 0,
            shiftFACE  = 1,
            shiftTRAY  = 2,
            shiftTOWER = 18,
            maskVIEW   = 1,
            maskFACE   = 1, 
            maskTRAY   = 0x0ff, 
            maskTOWER  = 0x0ff 
        };
    public:
        SortId(int tower, int tray, int face, int view) {
            m_sortId = (((maskTOWER & tower)<<shiftTOWER) | ((maskTRAY & tray)<<shiftTRAY) 
                | ((maskFACE & face)<< shiftFACE) | ((maskVIEW & view)<< shiftVIEW));
        }
        int getTower() { return ((m_sortId >> shiftTOWER) & maskTOWER);}
        int getTray()  { return ((m_sortId >> shiftTRAY) & maskTRAY); }
        int getFace()  { return ((m_sortId >> shiftFACE) & maskFACE); }
        int getView()  { return ((m_sortId >> shiftVIEW) & maskVIEW); }
        //operator int () const {return m_sortId; }
        bool operator < (const SortId& id) const { return (m_sortId < id.m_sortId); }

    private:
        int m_sortId;
    };


/**
* @class TkrTruncationInfo
*
* @brief  TkrTruncationInfo serves as a TDS container for a TkrTruncationMap.
*
* @author Leon Rochester
*
* $Header$
*/
    class TkrTruncationInfo : public DataObject {

    public:
        typedef std::map<SortId, TkrTruncatedPlane> TkrTruncationMap;

        /// Initializes the container
        TkrTruncationInfo() : m_nRCTrunc(0), m_nCCTrunc(0), m_done(false) {
            m_truncationMap.clear();
        }

        virtual ~TkrTruncationInfo() {
            m_truncationMap.clear();
        }

        /// Returns the TkrTruncationMap
        void setRCTrunc(int val) { m_nRCTrunc = val; }
        void setCCTrunc(int val) { m_nCCTrunc = val; }
        void setDone()           { m_done = true; }
        /// retrieves a pointer to the truncation map
        TkrTruncationMap* getTruncationMap() { return &m_truncationMap; }
        /// generates a key and adds the trucationItem to the map
        void addItem(int tower, int tray, int face, int view, TkrTruncatedPlane item) {
            SortId sortId(tower, tray, face, view);
            m_truncationMap[sortId] = item;
        }
        /// if true, there are truncated hits in the plane
        bool isTruncated() const { return ((m_nRCTrunc | m_nCCTrunc)!=0); }
        int  getNumRCTruncated() const { return m_nRCTrunc; }
        int  getNumCCTruncated() const { return m_nCCTrunc; }
        bool isDone()  const { return m_done; }

    private:

        /// map containing the truncation information
        TkrTruncationMap m_truncationMap;
        /// number of planes in the event with truncated read controller buffers
        int m_nRCTrunc;
        /// number of planes in the event with truncated cable controller buffers
        int m_nCCTrunc;
        /// if true, the truncation info for this event has been generated
        bool m_done;

    };
}

#endif
