## Process this file with automake to produce Makefile.in
#
# $Id$
#
# Author: Paul_Kunz@slac.stanford.edu
#

# The following set,  else would follow GNU conventions
AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = libgeomrep.la

libgeomrep_la_SOURCES = \
	ArcRep.cxx      \
	BoxRep.cxx      \
	ConeRep.cxx     \
	ConsRep.cxx     \
	HypeRep.cxx     \
	PointRep.cxx    \
	SpheRep.cxx     \
	TrackRep.cxx    \
	TubeRep.cxx     \
	TubsRep.cxx

libgeomrep_la_LDFLAGS = -version-info 1

INCLUDES = -I$(top_srcdir) \
	-I$(top_srcdir)/../CLHEP    \
	-I$(top_srcdir)/../geometry \
	-I$(top_srcdir)/../graphics \
	-I$(glast_prefix)/include
