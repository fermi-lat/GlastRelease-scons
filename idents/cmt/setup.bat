echo off
IF .%1==. (set CMTCONFIG=VisualC) ELSE set CMTCONFIG=%1
set CMTROOT=c:\packages\CMT\v1r5p1
set GLAST_SETTINGSROOT=c:\packages\external\glast_settings\v3
set GLAST_SETTINGSCONFIG=%CMTCONFIG%
set EXTLIBROOT=c:\packages\external\EXTLIB\v2r2
set EXTLIBCONFIG=%CMTCONFIG%
set IDENTSROOT=c:\packages\gaudi_dev\idents\v1r1
set IDENTSCONFIG=%CMTCONFIG%
set EXT_DIR=C:\extlib
set LHCXX_DIR=%EXT_DIR%
set ROOT_DIR=%EXT_DIR%\ROOT
set CLASSPATH=%CMTROOT%\java