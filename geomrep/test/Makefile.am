## Process this file with automake to produce Makefile.in
#
# $Id$
#
# Author: Paul_Kunz@slac.stanford.edu
#

# The following set,  else would follow GNU conventions
AUTOMAKE_OPTIONS = foreign

noinst_PROGRAMS = testgeom

testgeom_SOURCES = testgeom.cxx

INCLUDES = -I$(top_srcdir) \
	-I$(top_srcdir)/../geometry \
	-I$(top_srcdir)/../graphics \
	-I$(top_srcdir)/../CLHEP


test_LDADDS = \
	$(top_builddir)/src/libgeomrep.la

testgeom_LDADD = $(test_LDADDS)

TESTS = $(noinst_PROGRAMS)


