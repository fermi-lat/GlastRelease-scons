//     $Header$
//  author:  Toby Burnett


#ifndef WIN32

#include "Xostream.h"


#include <errno.h>
#include <strstream>
#include <iostream>
#include <cstring>
#include <cstdlib>

#include <cstdio>

//////////////////////////////////////////////////////////////////////
//	    private streambuf used by Xostream (from LAT)
class logbuf : public std::streambuf
{
    
public:
    enum { SIZE = 32 };

    logbuf (void);
    ~logbuf (void);

#if  !defined(__GNUG__)
    virtual int_type            overflow (int_type c /* = traits::eof()*/);
    virtual std::streamsize     xsputn (const char_type *string,
                                        std::streamsize length);
#else
    virtual int                 overflow (int c /* = EOF */);
    virtual int                 xsputn (const char *string, int length);
#endif
private:
    int                         m_pos;
#if  !defined(__GNUG__)
    char_type                   m_buffer [SIZE];
#else
    char                        m_buffer [SIZE];
#endif
};

//////////////////////////////////////////////////////////////////////
//	   Xostream implementations
//////////////////////////////////////////////////////////////////////
//		  constructor
Xostream::Xostream(const char* /*name*/)
: std::ostream(new logbuf)
{
}
//------------------------------------------------------------------------
//			  flush
void Xostream::flush()
{
    std::ostream::flush();
    // TODO: do I need this now?
}
//------------------------------------------------------------------------
//			     clear

void Xostream::clear()
{
    //TODO: figure out how to clear console
   flush();

}
//////////////////////////////////////////////////////////////////////////////
//	      logbuf implementations
//////////////////////////////////////////////////////////////////////////////

//<doc>         Initialize the log buffer.
logbuf::logbuf ()
    : m_pos (0)
{
    setp (0, 0);
    setg (0, 0, 0);
}

//<doc>         Flush the buffer.
logbuf::~logbuf (void)
{
#if  !defined(__GNUG__)
    overflow (traits_type::eof ());
#else
    overflow (EOF);
#endif
}

//<doc>         Handle buffer output.
//
//              Since the put area is empty we are in effect unbuffered, and
//              this method will be called every time there is output.  We
//              take the opportunity to adjust the line prefix here and then
//              buffer <<SIZE>> characters, or till end of the line, which
//              ever is less.  When the buffer must be flushed, it is sent
//              back to the <<logstream>>.
//
//<args>        ch                      Character to consume or end of file
//                                      indicator to flush. 
//
//<return>      End of file indicator on failure, zero on success.
#if  !defined(__GNUG__)
std::streambuf::int_type
logbuf::overflow (int_type ch)
#else
int
logbuf::overflow (int ch)
#endif
{
    if (
#if !defined(__GNUG__)
        ch != traits_type::eof ()
#else
        ch != EOF
#endif
        )
    {
        m_buffer [m_pos++] = (char) ch;

        // FIXME: should use the following, but MSVC++ does not have
        // `newline' in traits
        //   ch == traits_type::newline ()
        if (ch == '\n' || m_pos == SIZE)
        {
            std::cout.write (m_buffer, m_pos);
            m_pos = 0;
        }
    }
    else if (m_pos)
    {
        std::cout.write (m_buffer, m_pos);
        m_pos = 0;
    }

    return 0;
}

//<doc>         Handle lengthy buffer output.
//
//              This function permits faster output when a longer sequence of
//              known characters will be written.  Instead of writing every
//              character separately it prints out the whole thing.
//
//<args>        string                  The sequence of characters to output.
//              length                  The number of characters in the
//                                      sequence.
//
//<return>      The number of characters delivered to output; the value is
//              always equal to <<length>>.
#if !defined(__GNUG__)
std::streamsize
logbuf::xsputn (const char_type *string, std::streamsize length)
#else
int
logbuf::xsputn (const char *string, int length)
#endif
{
    for (int start = 0, end = 0; start < length; start = end) {
        // FIXME: should use: `&& string [end] != traits_type::newline ()'
        while (end < length && string [end] != '\n')
            ++end;

        if (end < length)
            // Include the newline
            ++end;

        // Copy either via buffer (if there was already something
        // there) or emit directly (if the buffer was empty or if
        // there are too many characters).  We use the buffer only on
        // the first string segment if the buffer is non-empty and on
        // the last segment if there isn't a trailing newline and the
        // segments fit into the buffer.
        if (m_pos) {
            // Use the buffer since there is already something there.
            while (m_pos < SIZE && start < end)
                m_buffer [m_pos++] = string [start++];
            if (m_pos == SIZE) {
                std::cout.write (m_buffer, m_pos);
                m_pos = 0;
            }
        }
        if (end - start < SIZE - m_pos)
            // Copy a trailing, not newline-terminated segment that is
            // short enough to fit to the buffer.  Note that the
            // previous `if' might have already flushed out a first
            // part of this string.
            while (start < end) m_buffer [m_pos++] = string [start++];
        else if (start < end)
            // Just write it out.  Either it is newline terminated or
            // it is the rest of a too long string, or it is part of a
            // multiline string.
            std::cout.write (string + start, end - start);
    }
    return length;
}
#endif
