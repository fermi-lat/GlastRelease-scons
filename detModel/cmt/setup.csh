# echo "Setting detModel v1r1 in /scratch/users/riccardo/glast-slac"

setenv CMTROOT /atl2lx/users/riccardo/CMT/v1r6
source ${CMTROOT}/mgr/setup.csh

set tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet setup -csh -pack=detModel -version=v1r1 -path=/scratch/users/riccardo/glast-slac $* >${tempfile}; source ${tempfile}; /bin/rm -f ${tempfile}
