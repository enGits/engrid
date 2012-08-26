@echo off
SETLOCAL EnableDelayedExpansion
rem =============== Developer configuration area ===============

rem VTK: environment variables for the chosen VTK installation
set VTKINCDIR=P:\ParaView-3.12.0-win64-fulldev\include\paraview-3.12
set VTKBINDIR=P:\ParaView-3.12.0-win64-fulldev\bin
set VTKLIBDIR=P:\ParaView-3.12.0-win64-fulldev\lib\paraview-3.12

rem Qt: batch file with environment variables for the chosen qt installation
set QTbatchfile=Q:\4.6.3_x64\qtvars.bat

rem NSIS: path to where NSIS is installed
set NSIS_DIR=P:\NSIS

rem BRLCAD: path to where BRLCAD-devel is installed
set BRLCAD_DIR=P:\BRLCADDEV\brlcad-7.22.0_dll_devel

rem ============= End Developer configuration area =============





rem ===================== Setup and copy =======================
rem go to where this script is located
%~d0
cd %~dp0
cd ..

rem go to the third_party64 folder 
IF NOT EXIST third_party64 mkdir third_party64
cd third_party64 || goto BADEND

rem Copy VTK installation
echo Copying VTK files, please wait...
IF NOT EXIST VTK mkdir VTK
cd VTK
IF NOT EXIST bin mkdir bin
IF NOT EXIST lib mkdir lib
IF NOT EXIST include mkdir include

xcopy /s "%VTKBINDIR%\*.dll" bin\ > NUL:
xcopy /s "%VTKLIBDIR%\*.*" lib\ > NUL:
xcopy /s "%VTKINCDIR%\*.*" include\ > NUL:

rem remove unwanted Qt*.dll files
del bin\Qt*.dll

rem going back to the third_party64 folder
cd ..


rem Copy Qt Installation
IF NOT EXIST Qt mkdir Qt
cd Qt

IF NOT EXIST bin mkdir bin
IF NOT EXIST lib mkdir lib
IF NOT EXIST src mkdir src
IF NOT EXIST include mkdir include
IF NOT EXIST plugins mkdir plugins
IF NOT EXIST mkspecs mkdir mkspecs
IF NOT EXIST tools mkdir tools
IF NOT EXIST translations mkdir translations

call %QTbatchfile%
echo Copying Qt files, please wait...

xcopy /s "%QTDIR%\bin\*.*" bin\ > NUL:
xcopy /s "%QTDIR%\lib\*.*" lib\ > NUL:
xcopy /s "%QTDIR%\src\*.*" src\ > NUL:
xcopy /s "%QTDIR%\include\*.*" include\ > NUL:
xcopy /s "%QTDIR%\plugins\*.*" plugins\ > NUL:
xcopy /s "%QTDIR%\mkspecs\*.*" mkspecs\ > NUL:
xcopy /s "%QTDIR%\tools\*.*" tools\ > NUL:
xcopy /s "%QTDIR%\translations\*.*" translations > NUL:

rem hack the qt.conf file
echo [Paths] > bin\qt.conf
echo Prefix = %CD% >> bin\qt.conf
echo Demos = demos >> bin\qt.conf
echo Examples = examples >> bin\qt.conf

rem create the qtvars.bat file
echo @echo off > bin\qtvars.bat
echo echo Setting QMAKESPEC to %QMAKESPEC% >> bin\qtvars.bat
echo set QMAKESPEC=%QMAKESPEC% >> bin\qtvars.bat
echo echo Setting QTDIR environment variable to %CD% >> bin\qtvars.bat
echo set QTDIR=%CD% >> bin\qtvars.bat
echo echo Putting Qt\bin in the current PATH environment variable. >> bin\qtvars.bat
echo set PATH=%CD%\bin;^%%PATH^%% >> bin\qtvars.bat
echo. >> bin\qtvars.bat

rem going back to the third_party64 folder
cd ..


rem Copy NSIS Installation
echo Copying NSIS files, please wait...
IF NOT EXIST NSIS mkdir NSIS
cd NSIS
xcopy /s "%NSIS_DIR%\*.*" . > NUL:

rem create the nsisvars.bat file
echo @echo off > nsisvars.bat
echo echo Putting NSIS in the current PATH environment variable. >> nsisvars.bat
echo set PATH=%CD%;^%%PATH^%% >> nsisvars.bat
echo. >> nsisvars.bat

rem going back to the third_party64 folder
cd ..


rem Copy BRLCAD Devel files
echo Copying BRLCAD Development files, please wait...
IF NOT EXIST BRLCAD mkdir BRLCAD
cd BRLCAD

IF NOT EXIST bin mkdir bin
IF NOT EXIST lib mkdir lib
IF NOT EXIST include mkdir include
IF NOT EXIST include\brlcad mkdir include\brlcad

xcopy /s "%BRLCAD_DIR%\bin64\*.*" bin\ > NUL:
xcopy /s "%BRLCAD_DIR%\lib64\*.*" lib\ > NUL:
xcopy /s "%BRLCAD_DIR%\include\*.*" include\brlcad\ > NUL:
xcopy /s "%BRLCAD_DIR%\COPYING" . > NUL:
xcopy /s "%BRLCAD_DIR%\readme.txt" . > NUL:
xcopy /s "%BRLCAD_DIR%\doc" . > NUL:

rem going back to the third_party folder
cd ..


rem =============== ALL DONE ================
:END
echo All done!
pause
exit 0

:BADEND
echo Unable to complete procedure.
pause
exit 1
