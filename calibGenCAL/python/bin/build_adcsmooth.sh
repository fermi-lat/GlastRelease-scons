#! /bin/bash
#$Header$
export -n DISPLAY


python ${CALIBGENCALROOT}/python/build_adcsmooth.py "$@"
