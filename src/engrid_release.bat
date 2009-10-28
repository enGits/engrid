REM Set up variables
set NSI_FILE=engrid.nsi

REM compile new (in release mode!)
REM C:\Qt\2009.03\bin\qtenv.bat
echo Setting up a MinGW/Qt only environment...
echo -- QTDIR set to C:\Qt\2009.03\qt
echo -- PATH set to C:\Qt\2009.03\qt\bin
echo -- Adding C:\Qt\2009.03\bin to PATH
echo -- Adding %SystemRoot%\System32 to PATH
echo -- QMAKESPEC set to win32-g++
set QTDIR=C:\Qt\2009.03\qt
set PATH=C:\Qt\2009.03\qt\bin
set PATH=%PATH%;C:\Qt\2009.03\bin;C:\Qt\2009.03\mingw\bin
set PATH=%PATH%;%SystemRoot%\System32
set QMAKESPEC=win32-g++

REM clean up
qmake
IF ERRORLEVEL 1 EXIT /B
C:/Qt/2009.03/mingw/bin/mingw32-make distclean
IF ERRORLEVEL 1 EXIT /B

REM build
qmake
IF ERRORLEVEL 1 EXIT /B
C:/Qt/2009.03/mingw/bin/mingw32-make -f Makefile.Release
IF ERRORLEVEL 1 EXIT /B

REM create setup for new
"C:\Program Files\NSIS\makensis.exe" "%NSI_FILE%"
IF ERRORLEVEL 1 EXIT /B
