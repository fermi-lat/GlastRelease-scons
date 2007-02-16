#! /bin/sh
#$Header$
export -n DISPLAY


python ${CALIBGENCALROOT}/python/mpdApplyCalibGain.py "$@"
