#!/bin/sh
mysqldump4 -h glastDB.slac.stanford.edu --result-file=staticDump.sql -p mood_test Parameter_class FSW_class Offline_class Ancillary_class FClass_to_PClass OClass_to_AClass PClass_to_AClass Instrument
