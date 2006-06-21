#! /bin/sh
#$Header$

PYTHONPATH=${CALIBGENCALROOT}/python/lib:${PYTHONPATH}
export PYTHONPATH

python ${CALIBGENCALROOT}/python/dumpROOTPlots.py "$@"

