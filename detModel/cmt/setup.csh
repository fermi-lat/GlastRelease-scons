# echo "Setting detModel v2r1 in /afs/slac.stanford.edu/u/ey/jrb/glast/jrbPack"

setenv CMTROOT /u/ey/jrb/glast/CMT/v1r6p1
source ${CMTROOT}/mgr/setup.csh

set tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet setup -csh -pack=detModel -version=v2r1 -path=/afs/slac.stanford.edu/u/ey/jrb/glast/jrbPack $* >${tempfile}; source ${tempfile}; /bin/rm -f ${tempfile}
