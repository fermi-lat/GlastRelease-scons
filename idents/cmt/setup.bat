echo off
IF .%1==. (set CMTCONFIG=VisualC) ELSE set CMTCONFIG=%1
set CMTROOT=d:\packages\CMT\v1r7
set GLASTPOLICYROOT=d:\packages\gismosys_v4\GlastPolicy\v1
set GLASTPOLICYCONFIG=Win32Debug
set GAUDIPOLICYROOT=d:\packages\gaudi_v7\GaudiPolicy\v4r3
set GAUDIPOLICYCONFIG=Win32Debug
set IDENTSROOT=d:\packages\pdr_dev\idents\v2
set IDENTSCONFIG=Win32Debug
set BINDIR=Win32Debug
set CLASSPATH=c:\Program Files\Java Analysis Studio\lib\jas.jar;c:\Program Files\Java Analysis Studio\lib\xml4j.jar;%CMTROOT%\java
set include=C:\Program Files\Microsoft Visual Studio\VC98\atl\include;C:\Program Files\Microsoft Visual Studio\VC98\include;c:\program files\devstudio\vc\include;c:\program files\devstudio\vc\atl\include;c:\program files\devstudio\vc\mfc\include;%include%
set lib=C:\Program Files\Microsoft Visual Studio\VC98\lib;c:\program files\devstudio\vc\lib;c:\program files\devstudio\vc\mfc\lib;%lib%