#! /bin/bash
#$Header$
export -n DISPLAY


python ${CALIBGENCALROOT}/python/pedDiff.py "$@"
