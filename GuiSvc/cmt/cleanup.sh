tempfile=/tmp/setup$$
${CMTROOT}/mgr/cmt -quiet cleanup -sh -pack=GuiSvc -version=v1r3 -path=/a/surrey10/g.glast_users/glground/tlindner/packages3/pdrApp_v2 $* >$tempfile; . $tempfile; /bin/rm -f $tempfile

