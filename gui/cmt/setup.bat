echo off
IF .%1==. (set CMTCONFIG=VisualC) ELSE set CMTCONFIG=%1
set CMTROOT=c:\packages\CMT\v1r5p1
set GLAST_SETTINGSROOT=c:\packages\external\glast_settings\v3r1
set GLAST_SETTINGSCONFIG=%CMTCONFIG%
set CLHEPROOT=c:\packages\external\CLHEP\v1r6
set CLHEPCONFIG=%CMTCONFIG%
set EXTLIBROOT=c:\packages\external\EXTLIB\v2r3
set EXTLIBCONFIG=%CMTCONFIG%
set GUIROOT=c:\packages\gui_dev\gui\v1
set GUICONFIG=%CMTCONFIG%
set EXT_DIR=%GLAST_EXT%
set LHCXX_DIR=%EXT_DIR%
set ROOT_DIR=%EXT_DIR%\ROOT
set CLHEP_DIR=%LHCXX_DIR%\CLHEP\1.6.0.0
set CLASSPATH=%CMTROOT%\java