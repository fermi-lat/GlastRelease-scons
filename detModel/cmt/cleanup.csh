set tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet cleanup -csh -pack=detModel -version=HEAD -path=/scratch/users/riccardo/glast-slac $* >$tempfile; source $tempfile; /bin/rm -f $tempfile

