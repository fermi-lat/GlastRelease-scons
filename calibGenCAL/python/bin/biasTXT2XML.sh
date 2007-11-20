#! /bin/bash
#$Header$
export -n DISPLAY


python ${CALIBGENCALROOT}/python/biasTXT2XML.py "$@"
