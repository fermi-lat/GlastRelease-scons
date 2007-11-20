#! /bin/bash
#$Header$
export -n DISPLAY


python ${CALIBGENCALROOT}/python/mpdApplyCalibGain.py "$@"
