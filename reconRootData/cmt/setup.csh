# echo "Setting reconRoot v0 in I:/packages/Root2IDL"

setenv CMTROOT I:\packages/CMT/v1r5p1
source ${CMTROOT}/mgr/setup.csh

set tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet setup -csh -pack=reconRoot -version=v0 -path=I:/packages/Root2IDL $* >$tempfile; source $tempfile; /bin/rm -f $tempfile
