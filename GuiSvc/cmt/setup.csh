# echo "Setting GuiSvc v1r3 in /a/surrey10/g.glast_users/glground/tlindner/packages3/pdrApp_v2"

setenv CMTROOT /afs/slac.stanford.edu/g/glast/applications/CMT/v1r6p1
source ${CMTROOT}/mgr/setup.csh

set tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet setup -csh -pack=GuiSvc -version=v1r3 -path=/a/surrey10/g.glast_users/glground/tlindner/packages3/pdrApp_v2 $* >${tempfile}; source ${tempfile}; /bin/rm -f ${tempfile}
