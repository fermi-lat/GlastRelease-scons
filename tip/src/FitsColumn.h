/** \file FitsColumn.h
    \brief FITS-specific interface to columnar data.
    \author James Peachey, HEASARC
*/
#ifndef tip_FitsColumn_h
#define tip_FitsColumn_h

#include <cassert>
#include <cstdlib>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>

#include "fitsio.h"

#include "FitsTable.h"
#include "FitsPrimProps.h"
#include "tip/IColumn.h"
#include "tip/TipException.h"
#include "tip/tip_types.h"

namespace tip {

  template <typename T>
  class FitsColumn : public IColumn {
    public:
      FitsColumn(FitsTable * ext, const std::string & id, FieldIndex_t field_index);

      virtual ~FitsColumn() throw() {}

      virtual void get(Index_t record_index, double & dest) const { getScalar(record_index, dest); }
      virtual void get(Index_t record_index, float & dest) const { getScalar(record_index, dest); }
      virtual void get(Index_t record_index, char & dest) const { getScalar(record_index, dest); }
      virtual void get(Index_t record_index, signed char & dest) const { getScalar(record_index, dest); }
      virtual void get(Index_t record_index, signed short & dest) const { getScalar(record_index, dest); }
      virtual void get(Index_t record_index, signed int & dest) const { getScalar(record_index, dest); }
      virtual void get(Index_t record_index, signed long & dest) const { getScalar(record_index, dest); }
      virtual void get(Index_t record_index, unsigned char & dest) const { getScalar(record_index, dest); }
      virtual void get(Index_t record_index, unsigned short & dest) const { getScalar(record_index, dest); }
      virtual void get(Index_t record_index, unsigned int & dest) const { getScalar(record_index, dest); }
      virtual void get(Index_t record_index, unsigned long & dest) const { getScalar(record_index, dest); }

      virtual void get(Index_t record_index, std::vector<double> & dest) const { getVector(record_index, dest); }
      virtual void get(Index_t record_index, std::vector<float> & dest) const { getVector(record_index, dest); }
      virtual void get(Index_t record_index, std::vector<char> & dest) const { getVector(record_index, dest); }
      virtual void get(Index_t record_index, std::vector<signed char> & dest) const { getVector(record_index, dest); }
      virtual void get(Index_t record_index, std::vector<signed short> & dest) const { getVector(record_index, dest); }
      virtual void get(Index_t record_index, std::vector<signed int> & dest) const { getVector(record_index, dest); }
      virtual void get(Index_t record_index, std::vector<signed long> & dest) const { getVector(record_index, dest); }
      virtual void get(Index_t record_index, std::vector<unsigned char> & dest) const { getVector(record_index, dest); }
      virtual void get(Index_t record_index, std::vector<unsigned short> & dest) const { getVector(record_index, dest); }
      virtual void get(Index_t record_index, std::vector<unsigned int> & dest) const { getVector(record_index, dest); }
      virtual void get(Index_t record_index, std::vector<unsigned long> & dest) const { getVector(record_index, dest); }

      virtual void set(Index_t record_index, const double & src) { setScalar(record_index, src); }
      virtual void set(Index_t record_index, const float & src) { setScalar(record_index, src); }
      virtual void set(Index_t record_index, const char & src) { setScalar(record_index, src); }
      virtual void set(Index_t record_index, const signed char & src) { setScalar(record_index, src); }
      virtual void set(Index_t record_index, const signed short & src) { setScalar(record_index, src); }
      virtual void set(Index_t record_index, const signed int & src) { setScalar(record_index, src); }
      virtual void set(Index_t record_index, const signed long & src) { setScalar(record_index, src); }
      virtual void set(Index_t record_index, const unsigned char & src) { setScalar(record_index, src); }
      virtual void set(Index_t record_index, const unsigned short & src) { setScalar(record_index, src); }
      virtual void set(Index_t record_index, const unsigned int & src) { setScalar(record_index, src); }
      virtual void set(Index_t record_index, const unsigned long & src) { setScalar(record_index, src); }

      virtual void set(Index_t record_index, const std::vector<double> & src) { setVector(record_index, src); }
      virtual void set(Index_t record_index, const std::vector<float> & src) { setVector(record_index, src); }
      virtual void set(Index_t record_index, const std::vector<char> & src) { setVector(record_index, src); }
      virtual void set(Index_t record_index, const std::vector<signed char> & src) { setVector(record_index, src); }
      virtual void set(Index_t record_index, const std::vector<signed short> & src) { setVector(record_index, src); }
      virtual void set(Index_t record_index, const std::vector<signed int> & src) { setVector(record_index, src); }
      virtual void set(Index_t record_index, const std::vector<signed long> & src) { setVector(record_index, src); }
      virtual void set(Index_t record_index, const std::vector<unsigned char> & src) { setVector(record_index, src); }
      virtual void set(Index_t record_index, const std::vector<unsigned short> & src) { setVector(record_index, src); }
      virtual void set(Index_t record_index, const std::vector<unsigned int> & src) { setVector(record_index, src); }
      virtual void set(Index_t record_index, const std::vector<unsigned long> & src) { setVector(record_index, src); }

      virtual void get(Index_t record_index, std::string & dest) const {
        char * buf = 0;
        try {
          // For strings, make a buffer to hold the value.
          buf = new char[m_width + 1];

          // Set the buffer contents to 0.
          memset(buf, '\0', m_width + 1);

          // Now call fitsio to fill the buffer with the string.
          getScalar(record_index, buf);

          // Copy this to the input.
          dest = buf;
        } catch (...) {
          delete [] buf;
          throw;
        }
        delete [] buf;
      }

      // TODO: Specialize this so that vectors of strings can be supported.
      // virtual void get(Index_t record_index, std::vector<std::string> & dest) const
      virtual void get(Index_t , std::vector<std::string> &) const
        { throw TipException("FitsColumn::get(Index_t, std::vector<std::string> &) not yet supported"); }

      virtual void set(Index_t record_index, const std::string & src) {
        setScalar(record_index, src.c_str());
      }

      // TODO: Specialize this so that vectors of strings can be supported.
      // virtual void set(Index_t record_index, const std::vector<std::string> & src)
      virtual void set(Index_t , const std::vector<std::string> &)
        { throw TipException("FitsColumn::set(Index_t, const std::vector<std::string> &) not yet supported"); }

      // Specializations for bool. This is done instead of specializing the templates to avoid complexity and
      // platform-specific code.
      virtual void get(Index_t record_index, bool & dest) const {
        if (!m_scalar) throw TipException("FitsColumn::get(Index_t, bool &) was called but field is not a scalar");
        int status = 0;
        char tmp_dest = 0;
        fits_read_col(m_ext->getFp(), FitsPrimProps<bool>::dataTypeCode(), m_field_index, record_index + 1, 1, m_repeat,
          0, &tmp_dest, 0, &status);
        if (0 != status) throw FitsTipException(status, "FitsColumn::get(Index_t, bool &) failed to read scalar cell value");
        dest = tmp_dest;
      }

      virtual void get(Index_t record_index, std::vector<bool> & dest) const {
        if (m_scalar) throw TipException("FitsColumn::get(Index_t, std::vector<bool> &) was called but field is not a vector");
        int status = 0;
        long num_els = getNumElements(record_index);
        char * tmp_dest = new char[num_els];
        fits_read_col(m_ext->getFp(), FitsPrimProps<bool>::dataTypeCode(), m_field_index, record_index + 1, 1, num_els,
          0, tmp_dest, 0, &status);
        if (0 != status) {
          delete [] tmp_dest;
          throw FitsTipException(status, "FitsColumn::get(Index_t, std::vector<bool> &) failed to read vector cell value");
        }
        dest.assign(tmp_dest, tmp_dest + num_els);
        delete [] tmp_dest;
      }

      virtual void set(Index_t record_index, const bool & src) {
        if (!m_scalar) throw TipException("FitsColumn::set(Index_t, const bool &) called but field is not a scalar");
        int status = 0;
        char tmp_src = src;
        fits_write_col(m_ext->getFp(), FitsPrimProps<bool>::dataTypeCode(), m_field_index, record_index + 1, 1, m_repeat,
          const_cast<void *>(static_cast<const void *>(&tmp_src)), &status);
        if (0 != status) throw FitsTipException(status, "FitsColumn::set(Index_t, const bool &) failed to write scalar cell value");
      }

      virtual void set(Index_t record_index, const std::vector<bool> & src) {
        if (m_scalar)
          throw TipException("FitsColumn::set(Index_t, const std::vector<bool> &) called but field is not a vector");
        int status = 0;
        long num_els = src.size();

        if (!m_var_length && num_els > m_repeat) {
          std::ostringstream os;
          os << "FitsColumn::set(Index_t, const std::vector<bool> &) attempted to write " << num_els <<
            " elements into a cell of size " << m_repeat;
          throw TipException(os.str());
        }

        char * tmp_src = new char[num_els];
        for (long ii = 0; ii < num_els; ++ii) tmp_src[ii] = src[ii];
        fits_write_col(m_ext->getFp(), FitsPrimProps<bool>::dataTypeCode(), m_field_index, record_index + 1, 1, num_els,
          tmp_src, &status);
        delete [] tmp_src;
        if (0 != status)
          throw FitsTipException(status, "FitsColumn::set(Index_t, const std::vector<bool> &) failed to write vector cell value");
      }

      /** \brief Copy a cell from another column to this column.
          \param src Pointer to the source column.
          \param src_index Index of the cell in the source column.
          \param dest_index Index of the cell in the destination column (this column).
      */
      virtual void copy(const IColumn * src, Index_t src_index, Index_t dest_index) {
        if (m_scalar) {
          T val;
          src->get(src_index, val);
          set(dest_index, val);
        } else {
          std::vector<T> val;
          src->get(src_index, val);
          set(dest_index, val);
        }
      }

      /** \brief Return a flag indicating whether this column holds scalar data.
      */
      virtual bool isScalar() const { return m_scalar; }

      /** \brief Return the name of the particular column implementation.
      */
      virtual const std::string implementation() const  { return "FITS"; }

      /** \brief Get number of elements in the given cell.
          \param record_index The record number identifying the cell.
      */
      virtual long getNumElements(Index_t record_index = 0) const {
        if (!m_var_length) return m_repeat;
        int status = 0;
        long num_els = 0;
        // Get number of elements in this particular field.
        fits_read_descript(m_ext->getFp(), m_field_index, record_index + 1, &num_els, 0, &status);
        if (0 != status) throw FitsTipException(status, "FitsColumn::getNumElements failed to get size of variable length cell");

        return num_els;
      }

      /** \brief Set number of elements in the given cell.
      */
      virtual void setNumElements(Index_t num_elements) {
        if (m_var_length) throw TipException("FitsColumn::setNumElements cannot change the width of variable length column");
        int status = 0;
        fits_modify_vector_len(m_ext->getFp(), m_field_index, num_elements, &status);
        if (0 != status) throw FitsTipException(status, "FitsColumn::setNumElements failed to modify field");

        // Update column information.
        m_repeat = num_elements;
        if (1 == m_repeat && !m_var_length) m_scalar = true; else m_scalar = false;
      }

    private:
      template <typename U>
      void getScalar(Index_t record_index, U & dest) const {
        // Prevent accidental calling for bool or string. The optimizer will swallow this.
        assert(typeid(U) != typeid(bool) && typeid(U) != typeid(std::string));
        if (!m_scalar) throw TipException("FitsColumn::getScalar was called but field is not a scalar");
        int status = 0;
        fits_read_col(m_ext->getFp(), FitsPrimProps<U>::dataTypeCode(), m_field_index, record_index + 1, 1, m_repeat,
          0, &dest, 0, &status);
        if (0 != status) throw FitsTipException(status, "FitsColumn::getScalar failed to read scalar cell value");
      }

      template <typename U>
      void getVector(Index_t record_index, std::vector<U> & dest) const {
        // Prevent accidental calling for bool or string. The optimizer will swallow this.
        assert(typeid(U) != typeid(bool) && typeid(U) != typeid(std::string));
        if (m_scalar) throw TipException("FitsColumn::getVector was called but field is not a vector");
        int status = 0;
        long num_els = getNumElements(record_index);
        dest.resize(num_els);
        U * dest_begin = &dest.front();
        fits_read_col(m_ext->getFp(), FitsPrimProps<U>::dataTypeCode(), m_field_index, record_index + 1, 1, num_els, 0,
          dest_begin, 0, &status);
        if (0 != status) throw FitsTipException(status, "FitsColumn::getVector failed to read vector cell value");
      }

      template <typename U>
      void setScalar(Index_t record_index, const U & dest) {
        // Prevent accidental calling for bool or string. The optimizer will swallow this.
        assert(typeid(U) != typeid(bool) && typeid(U) != typeid(std::string));
        if (!m_scalar) throw TipException("FitsColumn::setScalar called but field is not a scalar");
        int status = 0;
        if (m_ext->readOnly()) throw TipException("FitsColumn::setScalar called for a read-only file");
        fits_write_col(m_ext->getFp(), FitsPrimProps<U>::dataTypeCode(), m_field_index, record_index + 1, 1, m_repeat,
          const_cast<void *>(static_cast<const void *>(&dest)), &status);
        if (0 != status) throw FitsTipException(status, "FitsColumn::setScalar failed to write scalar cell value");
      }

      template <typename U>
      void setVector(Index_t record_index, const std::vector<U> & src) {
        // Prevent accidental calling for bool or string. The optimizer will swallow this.
        assert(typeid(U) != typeid(bool) && typeid(U) != typeid(std::string));
        if (m_scalar) throw TipException("FitsColumn::setVector called but field is not a vector");
        if (m_ext->readOnly()) throw TipException("FitsColumn::setVector called for a read-only file");
        int status = 0;
        long num_els = src.size();

        if (!m_var_length && num_els > m_repeat) {
          std::ostringstream os;
          os << "FitsColumn::setVector attempted to write " << num_els << " elements into a cell of size " << m_repeat;
          throw TipException(os.str());
        }

        const U * src_begin = &src.front();
        fits_write_col(m_ext->getFp(), FitsPrimProps<U>::dataTypeCode(), m_field_index, record_index + 1, 1, num_els,
          const_cast<void *>(static_cast<const void *>(src_begin)), &status);
        if (0 != status) throw FitsTipException(status, "FitsColumn::setVector failed to write vector cell value");
      }

      FitsTable * m_ext;
      FieldIndex_t m_field_index;
      long m_repeat;
      long m_width;
      int m_type_code;
      bool m_var_length;
      bool m_scalar;
  };

  template <typename T>
  inline FitsColumn<T>::FitsColumn(FitsTable * ext, const std::string & id, FieldIndex_t field_index): IColumn(id),
    m_ext(ext), m_field_index(field_index), m_repeat(0), m_width(0), m_type_code(0), m_var_length(false), m_scalar(false) {
    // Determine characteristics of this column.
    int status = 0;
    fits_get_coltype(m_ext->getFp(), m_field_index, &m_type_code, &m_repeat, &m_width, &status);
    if (0 != status) throw FitsTipException(status, "FitsColumn::FitsColumn failed to get information about field");

    // Handle special case of strings, for which the info returned by fits_get_coltype means something different.
    if (TSTRING == m_type_code) {
      m_repeat /= m_width;
      if (0 == m_repeat) m_repeat = 1;
    }

    // Handle variable length columns.
    if (m_type_code < 0) {
      m_type_code *= -1;
      m_var_length = true;
    }

    // Determine whether this is a scalar column.
    if (1 == m_repeat && !m_var_length) m_scalar = true;

    // Get units.
    std::ostringstream os;
    os << "TUNIT" << m_field_index;
    char units[FLEN_CARD] = "";
    fits_read_key(m_ext->getFp(), TSTRING, const_cast<char *>(os.str().c_str()), units, 0, &status);
    if (0 == status) m_units = units;
    else if (KEY_NO_EXIST != status) throw FitsTipException(status, "FitsColumn::FitsColumn failed to get units of field");
  }

}

#endif
