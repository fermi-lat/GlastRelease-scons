# echo "Setting detModel v2r1 in /afs/slac.stanford.edu/u/ey/jrb/glast/jrbPack"

CMTROOT=/u/ey/jrb/glast/CMT/v1r6p1; export CMTROOT
. ${CMTROOT}/mgr/setup.sh

tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet setup -sh -pack=detModel -version=v2r1 -path=/afs/slac.stanford.edu/u/ey/jrb/glast/jrbPack $* >${tempfile}; . ${tempfile}; /bin/rm -f ${tempfile}
