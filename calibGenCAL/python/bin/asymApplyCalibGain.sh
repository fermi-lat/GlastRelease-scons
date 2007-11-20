#! /bin/bash
#$Header$
export -n DISPLAY


python ${CALIBGENCALROOT}/python/asymApplyCalibGain.py "$@"
