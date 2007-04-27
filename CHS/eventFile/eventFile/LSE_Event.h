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

  /**
   * @brief Represents an event read in from a .evt file
   *
   * This collection of event information is produced only by an LSEReader and is
   * not to be modified or copied outside of that class.  Client software will refer
   * to instances via the smart-pointer type EventPtr which does reference counting.
   */
  class LSE_Event: private boost::noncopyable {

    // data members read from the .evt file
    const LSE_Context m_ctx;
    const EBF_Data    m_ebf;
    const LSE_Info::InfoType m_itype;
    const LPA_Info    m_pinfo;
    const LCI_ACD_Info m_ainfo;
    const LCI_CAL_Info m_cinfo;
    const LCI_TKR_Info m_tinfo;
    const LSE_Keys::KeysType m_ktype;
    const LPA_Keys     m_pakeys;
    const LCI_Keys     m_cikeys;

    // only an LSEReader is allowed to make instances of us.
    friend class LSEReader;
    LSE_Event( const LSE_Context& c, const EBF_Data& e, LSE_Info::InfoType i, const LPA_Info& pi,
	       const LCI_ACD_Info& ai, const LCI_CAL_Info& ci, const LCI_TKR_Info& ti,
	       LSE_Keys::KeysType k, const LPA_Keys& pk, const LCI_Keys& ck )
      : m_ctx( c ), m_ebf( e ), m_itype( i ), m_pinfo( pi ), m_ainfo( ai ), m_cinfo( ci ), 
	m_tinfo( ti ), m_ktype( k ), m_pakeys( pk ), m_cikeys( ck ) {}

  public:

    // need a public dtor
    ~LSE_Event() {}

    // accessors for the data members
    const LSE_Context& ctx() const { return m_ctx; }
    const EBF_Data&    ebf() const { return m_ebf; }
    LSE_Info::InfoType infotype() const { return m_itype; }
    LSE_Keys::KeysType keystype() const { return m_ktype; }
    const LPA_Info* pinfo() const { return ( m_itype == LSE_Info::LPA ) ? &m_pinfo : 0; }
    const LCI_ACD_Info* ainfo() const { return ( m_itype == LSE_Info::LCI_ACD ) ? &m_ainfo : 0; }
    const LCI_CAL_Info* cinfo() const { return ( m_itype == LSE_Info::LCI_CAL ) ? &m_cinfo : 0; }
    const LCI_TKR_Info* tinfo() const { return ( m_itype == LSE_Info::LCI_TKR ) ? &m_tinfo : 0; }
    const LPA_Keys* pakeys() const { return ( m_ktype == LSE_Keys::LPA ) ? &m_pakeys : 0; }
    const LCI_Keys* cikeys() const { return ( m_ktype == LSE_Keys::LCI ) ? &m_cikeys : 0; }
    
  };

};

/// everybody except LSEReader will use this smart pointer type to refer to LSE_Event's
typedef boost::shared_ptr< const eventFile::LSE_Event > EventPtr;

#endif // EVENTFILE_LSE_EVENT_H
