# echo "Setting DetDisp v1 in c:/packages/gui_dev"

setenv CMTROOT c:\packages/CMT/v1r5p1
source ${CMTROOT}/mgr/setup.csh

set tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet setup -csh -pack=DetDisp -version=v1 -path=c:/packages/gui_dev $* >$tempfile; source $tempfile; /bin/rm -f $tempfile
