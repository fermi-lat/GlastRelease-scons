#ifndef __CHECKSUM_H__rLCZxztLIzJGvRZ5GvXrYQWZyf7VMjyy7LDu
#define __CHECKSUM_H__rLCZxztLIzJGvRZ5GvXrYQWZyf7VMjyy7LDu

/**
 * @class checkSum
 *
 * @brief  checkSum is a simple class that plugs itself on RootTupleSvc to
 * calculate a checksum of the bytes being written to the MeritTuple ntuple.
 *
 * @author Michael Kuss
 *
 * $Header$
*/
 
#include "TTree.h"

#include <fstream>
#include <string>
#include <vector>


class checkSum {
 public:
    /// constructor taking the filename to store the checksum
    checkSum(const std::string);

    ~checkSum() {
        if ( m_is_open )
            m_out.close();
    }

    /// computes and writes the checksum to the file
    void write(TTree*);

    bool is_open() { return m_is_open; }
    bool bad() { return m_bad; }

private:

    /// computes a simple checksum
    unsigned long simple(std::vector<unsigned char>*);

    std::string m_fileName;

    /// ofstream of m_fileName
    std::ofstream m_out;

    /// set this flag if m_out is open
    bool m_is_open;

    /// set this flag if the file cannot be opened
    bool m_bad;
};

#endif
