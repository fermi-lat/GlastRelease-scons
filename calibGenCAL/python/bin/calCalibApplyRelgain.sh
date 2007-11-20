#! /bin/bash
#$Header$
export -n DISPLAY


python ${CALIBGENCALROOT}/python/calCalibApplyRelgain.py "$@"
