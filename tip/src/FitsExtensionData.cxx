/** \file FitsExtensionData.cxx

    \brief FITS-specific low level abstract interface to all extension types.
    \author James Peachey, HEASARC
*/

#include "FitsExtensionData.h"
#include "FitsExtensionUtils.h"
#include "FitsTabularData.h"
#include "table/IData.h"
#include "table/IHeaderData.h"

namespace table {

  // Construct without opening the file.
  FitsExtensionData::FitsExtensionData(const std::string & file_name, const std::string & ext_name):
    m_fits_utils(0), m_header(0), m_data(0) {
    try {
      m_fits_utils = new FitsExtensionUtils(file_name, ext_name);
      m_header = m_fits_utils->getHeaderData();
      m_data = m_fits_utils->getTabularData();
    } catch(...) {
//      delete m_data; m_data = 0;
//      delete m_header; m_header = 0;
      delete m_fits_utils; m_fits_utils = 0;
      throw;
    }
  }

  // Close file automatically while destructing.
  FitsExtensionData::~FitsExtensionData() { delete m_data; delete m_header; delete m_fits_utils; }

  IHeaderData * FitsExtensionData::getHeaderData() { return m_header; }

  ITabularData * FitsExtensionData::getTabularData() { return dynamic_cast<FitsTabularData *>(m_data); }

}
