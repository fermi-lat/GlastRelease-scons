# echo "Setting reconRootData v1r1 in /a/surrey10/g.glast_users/glground/tlindner/packages3/glastgaudi_test"

CMTROOT=/afs/slac.stanford.edu/g/glast/applications/CMT/v1r6p1; export CMTROOT
. ${CMTROOT}/mgr/setup.sh

tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet setup -sh -pack=reconRootData -version=v1r1 -path=/a/surrey10/g.glast_users/glground/tlindner/packages3/glastgaudi_test $* >${tempfile}; . ${tempfile}; /bin/rm -f ${tempfile}
