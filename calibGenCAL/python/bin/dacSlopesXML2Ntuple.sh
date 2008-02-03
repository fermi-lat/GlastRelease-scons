#! /bin/bash
#$Header$
export -n DISPLAY

python ${CALIBGENCALROOT}/python/dacSlopesXML2Ntuple.py "$@"

