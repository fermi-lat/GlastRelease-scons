#! /bin/bash
#$Header$
export -n DISPLAY


python ${CALIBGENCALROOT}/python/diff_LATC_Cal_Thresh.py "$@"
