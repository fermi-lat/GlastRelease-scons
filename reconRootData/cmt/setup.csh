# echo "Setting reconRootData v1r1 in /a/surrey10/g.glast_users/glground/tlindner/packages3/glastgaudi_test"

setenv CMTROOT /afs/slac.stanford.edu/g/glast/applications/CMT/v1r6p1
source ${CMTROOT}/mgr/setup.csh

set tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet setup -csh -pack=reconRootData -version=v1r1 -path=/a/surrey10/g.glast_users/glground/tlindner/packages3/glastgaudi_test $* >${tempfile}; source ${tempfile}; /bin/rm -f ${tempfile}
