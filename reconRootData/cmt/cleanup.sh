tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet cleanup -sh -pack=reconRootData -version=v1r1 -path=/a/surrey10/g.glast_users/glground/tlindner/packages3/glastgaudi_test $* >$tempfile; . $tempfile; /bin/rm -f $tempfile

