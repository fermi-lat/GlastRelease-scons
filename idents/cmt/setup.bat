echo off
IF .%1==. (set CMTCONFIG=VisualC) ELSE set CMTCONFIG=%1
set CMTROOT=c:\packages\CMT\v1r5p1
set GLAST_SETTINGSROOT=c:\packages\external\glast_settings\v2
set GLAST_SETTINGSCONFIG=%CMTCONFIG%
set EXTLIBROOT=c:\packages\external\EXTLIB\v2r3
set EXTLIBCONFIG=%CMTCONFIG%
set IDENTSROOT=c:\packages\gaudi_dev\idents\v1
set IDENTSCONFIG=%CMTCONFIG%
set EXT_DIR=%GLAST_EXT%
set LHCXX_DIR=%EXT_DIR%
set ROOT_DIR=%EXT_DIR%\ROOT
set CLASSPATH=%CMTROOT%\java