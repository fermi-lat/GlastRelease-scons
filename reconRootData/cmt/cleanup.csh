set tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet cleanup -csh -pack=reconRootData -version=v1r1 -path=/a/surrey10/g.glast_users/glground/tlindner/packages3/glastgaudi_test $* >$tempfile; source $tempfile; /bin/rm -f $tempfile

