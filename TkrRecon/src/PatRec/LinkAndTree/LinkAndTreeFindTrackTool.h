/**
 * @class LinkAndTreeFindTrackTool
 *
 * @brief Implements a Gaudi Tool for find track candidates. This particular tool uses the 
 *        "Combo" method.
 *
 * @author The Tracking Software Group
 *
 * $Header$
 */

#ifndef LINKANDTREEFINDTRACKTOOL_H
#define LINKANDTREEFINDTRACKTOOL_H

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "TkrRecon/PatRec/ITkrFindTrackTool.h"
#include "TkrRecon/Services/TkrGeometrySvc.h"

class LinkAndTreeFindTrackTool : public AlgTool, virtual public ITkrFindTrackTool
{
public:
    /// Standard Gaudi Tool interface constructor
    LinkAndTreeFindTrackTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~LinkAndTreeFindTrackTool() {}

    /// @brief Method to find candidate tracks. Will retrieve the necessary information from
    ///        the TDS, including calorimeter energy, and then use TkrLinkAndTree to find all
    ///        possible track candidates. The resulting track candidate collection is then 
    ///        stored in the TDS for the next stage.
    StatusCode findTracks();

private:
    /// Pointer to the local Tracker geometry service
    TkrGeometrySvc* m_tkrGeo;

    /// Pointer to the Gaudi data provider service (interface to the TDS)
    DataSvc*        m_dataSvc;
};

#endif
