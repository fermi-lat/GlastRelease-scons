# echo "Setting GraphicsSvc v1 in d:/packages/gaudi_dev"

setenv CMTROOT d:/packages/CMT/v1r5p1
source ${CMTROOT}/mgr/setup.csh

set tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet setup -csh -pack=GraphicsSvc -version=v1 -path=d:/packages/gaudi_dev $* >$tempfile; source $tempfile; /bin/rm -f $tempfile
