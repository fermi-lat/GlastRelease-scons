#! /bin/sh
#$Header$
export -n DISPLAY


python ${CALIBGENCALROOT}/python/get_slac_calibdac.py "$@"
