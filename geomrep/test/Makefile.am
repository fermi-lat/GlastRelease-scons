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

# X11 options
# the following should be added to files including X headers
X_CFLAGS = @X_CFLAGS@

# the path and X libraries needed
X_LIBS   = @X_LIBS@ -lXm -lXt -lX11
X_EXTRA_LIBS = @X_EXTRA_LIBS@

INCLUDES = -I$(top_srcdir) \
	-I$(top_srcdir)/../geometry \
	-I$(top_srcdir)/../graphics \
	-I$(top_srcdir)/../CLHEP \
	-I$(glast_prefix)/include

test_LDADDS =						\
	$(top_builddir)/libgeomrep.la			\
	$(top_builddir)/../graphics/libmotif.la		\
	$(top_builddir)/../graphics/libgraphics.la	\
	$(top_builddir)/../geometry/libgeometry.la	\
	$(top_builddir)/../CLHEP/CLHEP/libCLHEP.a	\
	$(X_LIBS)

testgeom_LDADD = $(test_LDADDS)

TESTS = $(noinst_PROGRAMS)


