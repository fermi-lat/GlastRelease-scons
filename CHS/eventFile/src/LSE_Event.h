// -*- mode: c++ -*-
/** @file LSE_Event.h
 *  @brief Defines class LSE_Event and the typedef EventPtr
 */

#ifndef EVENTFILE_LSE_EVENT_H
#define EVENTFILE_LSE_EVENT_H

#include <boost/shared_ptr.hpp>
#include <boost/utility.hpp>
#include <boost/weak_ptr.hpp>

#include "eventFile/LSE_Context.h"
#include "eventFile/LSE_Info.h"
#include "eventFile/EBF_Data.h"
#include "eventFile/LSE_Keys.h"

namespace eventFile {
  class LSEReader;
}

/**
 * @brief Represents an event read in from a .evt file
 *
 * This collection of event information is produced only by an LSEReader and is
 * not to be modified or copied outside of that class.  Client software will refer
 * to instances via the smart-pointer type EventPtr which does reference counting.
 */
class LSE_Event: private boost::noncopyable {

  // data members read from the .evt file
  const eventFile::LSE_Context m_ctx;
  const eventFile::EBF_Data    m_ebf;
  const eventFile::LSE_Info::InfoType m_itype;
  const eventFile::LPA_Info    m_pinfo;
  const eventFile::LCI_ACD_Info m_ainfo;
  const eventFile::LCI_CAL_Info m_cinfo;
  const eventFile::LCI_TKR_Info m_tinfo;
  const eventFile::LSE_Keys::KeysType m_ktype;
  const eventFile::LPA_Keys     m_pakeys;
  const eventFile::LCI_Keys     m_cikeys;

public:
  LSE_Event( const eventFile::LSE_Context& c, 
	     const eventFile::EBF_Data& e, 
	     eventFile::LSE_Info::InfoType i, 
	     const eventFile::LPA_Info& pi,
	     const eventFile::LCI_ACD_Info& ai, 
	     const eventFile::LCI_CAL_Info& ci, 
	     const eventFile::LCI_TKR_Info& ti,
	     eventFile::LSE_Keys::KeysType k, 
	     const eventFile::LPA_Keys& pk, 
	     const eventFile::LCI_Keys& ck )
    : m_ctx( c ), m_ebf( e ), m_itype( i ), m_pinfo( pi ), m_ainfo( ai ), m_cinfo( ci ), 
      m_tinfo( ti ), m_ktype( k ), m_pakeys( pk ), m_cikeys( ck ) {}

  ~LSE_Event() {}

  // accessors for the data members
  eventFile::LSE_Context ctx() const { return m_ctx; }
  const eventFile::EBF_Data&    ebf() const { return m_ebf; }
  eventFile::LSE_Info::InfoType infotype() const { return m_itype; }
  eventFile::LSE_Keys::KeysType keystype() const { return m_ktype; }
  const eventFile::LPA_Info* pinfo() const { return ( m_itype == eventFile::LSE_Info::LPA ) ? &m_pinfo : 0; }
  const eventFile::LCI_ACD_Info* ainfo() const { return ( m_itype == eventFile::LSE_Info::LCI_ACD ) ? &m_ainfo : 0; }
  const eventFile::LCI_CAL_Info* cinfo() const { return ( m_itype == eventFile::LSE_Info::LCI_CAL ) ? &m_cinfo : 0; }
  const eventFile::LCI_TKR_Info* tinfo() const { return ( m_itype == eventFile::LSE_Info::LCI_TKR ) ? &m_tinfo : 0; }
  const eventFile::LPA_Keys* pakeys() const { return ( m_ktype == eventFile::LSE_Keys::LPA ) ? &m_pakeys : 0; }
  const eventFile::LCI_Keys* cikeys() const { return ( m_ktype == eventFile::LSE_Keys::LCI ) ? &m_cikeys : 0; }
    
};

/// all clients will refer to this class through the smart pointer
typedef boost::shared_ptr<LSE_Event> EventPtr;

#endif // EVENTFILE_LSE_EVENT_H
