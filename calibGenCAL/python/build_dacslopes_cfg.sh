#! /bin/sh
#$Header$

PYTHONPATH=${CALIBGENCALROOT}/python/lib:${PYTHONPATH}
export PYTHONPATH

python ${CALIBGENCALROOT}/python/build_dacslopes_cfg.py "$@"
