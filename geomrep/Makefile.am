## Process this file with automake to produce Makefile.in
#
# $Id$
#
# Author: Paul_Kunz@slac.stanford.edu
#

# The following set,  else would follow GNU conventions
AUTOMAKE_OPTIONS = foreign

EXTRA_DIST = config

SUBDIRS = geomrep src . test

lib_LTLIBRARIES = libgeomrep.la

libgeomrep_la_SOURCES = version.cxx

libgeomrep_la_LIBADD = src/libgeomrep.la

libgeomrep_la_LDFLAGS = -version-info 1

