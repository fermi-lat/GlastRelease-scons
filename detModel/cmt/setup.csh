# echo "Setting detModel HEAD in /scratch/users/riccardo/glast-slac"

setenv CMTROOT /atl2lx/users/riccardo/CMT/v1r6
source ${CMTROOT}/mgr/setup.csh

set tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet setup -csh -pack=detModel -version=HEAD -path=/scratch/users/riccardo/glast-slac $* >${tempfile}; source ${tempfile}; /bin/rm -f ${tempfile}
