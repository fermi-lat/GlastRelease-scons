# echo "Setting GuiSvc v1r3 in /a/surrey10/g.glast_users/glground/tlindner/packages3/pdrApp_v2"

CMTROOT=/afs/slac.stanford.edu/g/glast/applications/CMT/v1r6p1; export CMTROOT
. ${CMTROOT}/mgr/setup.sh

tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet setup -sh -pack=GuiSvc -version=v1r3 -path=/a/surrey10/g.glast_users/glground/tlindner/packages3/pdrApp_v2 $* >${tempfile}; . ${tempfile}; /bin/rm -f ${tempfile}
